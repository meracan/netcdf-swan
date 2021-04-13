import os
import json
import numpy as np
from tqdm import tqdm
from scipy.io import loadmat
from s3netcdf import NetCDF2D 
from datetime import datetime
import re

class NetCDFSWAN(NetCDF2D):
  """
    Creates and stores partitioned netcdf (.nc) files from SWAN data into a local cache and
    uploads them to an Amazon S3 bucket, using the s3-netcdf package.
    
    This is the Swan output structure. 
    The year and month folders are saved as keys.
    `-- SWAN_DATA
        |
        |-- 2005
        |-- 2006
        |-- 2007
        :   |-- 01
            |-- 02
            |-- 03
            :   `-- results
                    |-- HS.mat
                    |-- WIND.mat
                    |-- line_n.spc
                    |-- QP.mat
                    :

                    :
                    |-- hotspots.spc
                    `-- TPS.mat
            :
            |-- 11
            `-- 12
        :
        |-- 2016
        |-- 2017
        |
        `-- Mesh
            |-- .ele
            |-- .bot
            `-- .node

    The Mesh folder holds latitudes, longitudes, and bathymetry of all the nodes, and other triangle mesh information.
    The Matlab (.mat) files in the 'results' folder of each month holds variable information for all of the nodes,
    as well as spectra information.
    
    Parameters
    ----------
    swanFolder:str
      Path of SWAN_DATA
    obj:obj (see test folder)
      object:
      name : str, Name of cdffile
      cacheLocation: str,path
      bucket: Name of S3 bucket
      localOnly:bool
        To include or ignore s3 storage
        Default:True
      ncSize:float,optional
        Max partition size (MB)   
      nca:object
        metadata:obj
        dimensions:[obj]
          name:str
          value:int
        groups:[obj]
          name:str,
          variables:[obj]
            name:str
            type:str
            dimensions:[str],
            etc..
    logger:logging.getLogger()
      Log basic printout to log file
  """
  def __init__(self,obj,logger=None):
    super().__init__(obj)
    
    self.swanFolder=obj.get("swanFolder")
    
    self.logger     = logger
  
    # Initialize some of the info from nca to self.
    # If the info is missing, it will raise an exception.
    info        = self.info()
    self.ntime  = ntime = info['dimensions'].get('ntime')
    self.nnode  = info['dimensions'].get('nnode')
    self.nsnode = info['dimensions'].get('nsnode')
    self.npart = info['dimensions'].get('npart')
    
    # Initialize datetime array for the entire simulation
    startDate = info['metadata'].get('startDate')
    timeStep  = info['metadata'].get('timeStep(h)')
    self.startDate = startDate=np.datetime64(startDate)
    self.datetime  = startDate+np.arange(ntime)*np.timedelta64(timeStep, 'h')
    
    # Load station index table
    # Each station (.spc) file needs a specific Id
    # The ids are saved under BSCWANv5.stations.json
    stations=info['metadata'].get('stations')
    self.stations = stations    
    
    # Load variables for matlab output
    # This contains an array of all variables and only used to save the temporal axis data
    variables=info['metadata'].get('mvariables')
    self.variables = variables
    
    # Get the matlab variable name
    mtname={}
    for key in self.variables:
      variable=self.variables.get(key)
      matname=variable.get("matfile name")
      mtname[matname]=key
    self.mtname=mtname
    
    # Combine matfile names and variables of partition files
    pt_gvar = set(info['groups']['pt']['variables'].keys())
    pt_intr = list(set(dict(variables.items()).keys()).intersection(pt_gvar))
    self.ptList = [v["fileKey"] for k,v in variables.items() if k in pt_gvar] + pt_intr

    # Collect corrupted files. Read from json first, if available
    errorlistPath = os.path.join(self.folder, "error_list.json")
    if not os.path.isfile(errorlistPath):json.dump([], open(errorlistPath, "w+"))
    self.errorlist = json.load(open(errorlistPath))

    
    # Show progress bar
    self.showProgress=showProgress=obj.get("showProgress",False)
    self.pbar0=None 
    if showProgress:
      try:
        from tqdm import tqdm
        self.pbar0 = tqdm(total=1, position=1)
        self.pbar = tqdm(total=1, position=0)
      except Exception as err:
        import warnings
        warnings.warn("tqdm does not exist")
    
    # Avoid certain files
    self.blacklist=info['metadata'].get('blacklist',[])
    
    # Number of nodes to upload at a time
    self.gnode = obj.get("gnode",1000)
 
  @property
  def matFiles(self):
    """ Get .mat files only
    """
    swanFolder=self.swanFolder
    if swanFolder is None:raise Exception("swanFolder was not specified")
    files = NetCDFSWAN.getFiles(swanFolder)
    return list(filter(lambda file:file['ext']=='.mat' and file['name'] not in self.blacklist,files))
  
  @property
  def spcFiles(self):
    """ Get .spc files only
    """
    swanFolder=self.swanFolder
    if swanFolder is None:raise Exception("swanFolder was not specified")
    files = NetCDFSWAN.getFiles(swanFolder)
    files =list(filter(lambda file:file['ext']=='.spc',files))
    return  list(sorted(files,key=lambda file:(file['year'],file['month'],file['name'])))
 
  @staticmethod
  def prepareInputJSON(file,swanFolder, **kwargs):
    """ Prepare input json object to load into NetCDF2D.
        This merge a few json together and extract certain info into metadata.
    """
    obj=NetCDFSWAN.load(file)
    obj['swanFolder']=swanFolder
    jsonFolder=os.path.dirname(file)
    
    # Load and merge json files BCWABv5.variables.json/BCWABv5.specvariables.json
    obj["nca"]['groups']['s']['variables']=f1=NetCDFSWAN.load(os.path.join(jsonFolder,obj["nca"]['groups']['s']['variables']))
    obj["nca"]['groups']['t']['variables']=f2=NetCDFSWAN.load(os.path.join(jsonFolder,obj["nca"]['groups']['t']['variables']))
    obj["nca"]['groups']['spc']['variables']=f3=NetCDFSWAN.load(os.path.join(jsonFolder,obj["nca"]['groups']['spc']['variables']))
    obj["nca"]['groups']['pt']['variables']=f4=NetCDFSWAN.load(os.path.join(jsonFolder,obj["nca"]['groups']['pt']['variables']))
    
    # Keep variables since we need to extract "matlab name" keys
    variables={**f1,**f2}
    
    # Store "matlab name" into a dictionnary
    temp_={**variables,**f3} # Need to keep the spectra name
    temp={**temp_, **f4} # partition variable names

    obj['nca']['metadata']['mvariables']=temp
    
    # Get spectral stations metadata that contains the station name and id
    meta=NetCDFSWAN.getSpectralStationMetadata(swanFolder, **kwargs)

    obj['nca']['metadata']['stations']=meta['stations']
    
    if obj['nca']['dimensions']['nsnode']!=meta['nsnode']:
      raise Exception("Please check nsnode in json file. {} to {}".format(obj['nca']['dimensions']['nsnode'],meta['nsnode']))
    
    if obj['nca']['dimensions']['nstation']!=meta['nstation']:
      raise Exception("Please check nstation in json file. {} to {}".format(obj['nca']['dimensions']['nstation'],meta['nstation']))      
    
    return obj
    
  @staticmethod
  def printMatKeys(swanFolder,year=2014):
    """ Simple function to print the keys within a matlab file
    """
    files=NetCDFSWAN.getFiles(swanFolder)
    files=list(filter(lambda file:"year" in file and file["year"]==year and "month" in file and file["month"]==1 and file['ext']=='.mat',files))
    for file in files:
      print(NetCDFSWAN.load(file["path"]).keys())
  
  @staticmethod
  def printSpcShape(swanFolder,year=2014):
    """ Simple function to print shape of the spectra files
    """    
    files=NetCDFSWAN.getFiles(swanFolder)
    files=list(filter(lambda file:"year" in file and file["year"]==year and "month" in file and file["month"]==1 and file['ext']=='.spc',files))
    for file in files:
      print("{} - {}".format( file["name"] , NetCDFSWAN.load(file["path"])['spectra'].shape))    
    
  
  @staticmethod
  def getSpectralStationMetadata(swanFolder,print_meta=False,year=2014,month=1):
    """ Extract spectralStation MetaData
        Note that this is project specific. 
        The code is using one of the folder output to extract the .spc files.
        The code calculates the number of .spc files/stations and the number of points within the file.
    """
    files=NetCDFSWAN.getFiles(swanFolder)
    
    files=list(filter(lambda file:"year"in file and "month" in file and file["year"]==year and file["month"]==month and file['ext']=='.spc',files))
    
    files=sorted(files,key=lambda x:x["name"])
    stations={}
    isnode=0
    for i,file in enumerate(files):
      nlatlng=NetCDFSWAN.load(file['path'],return_metadata=True)['nlatlng']    
      stations[file['name']]={"start":isnode,"end":isnode+nlatlng,"id":i}
      isnode += nlatlng
    
    nsnode=isnode
    if print_meta:
      print("Number of stations:{}".format(len(files)))
      print("Number of station points: {}".format(nsnode))
    
    return {"stations":stations,"nstation":len(files),"nsnode":nsnode}
    
  @staticmethod
  def getFiles(swanFolder):
    """ Get all output files from swanFolder and save it in a dictionnary.
        
        Returns
        -------
        list:object
          {"year":int,"month":int,"path":str,"name":str,"ext":str}
        
        Example:
          {"year":2004,"month":1,"path":"2004/01/results/HS.mat","name":"HS","ext":".mat"}
          {"path":"Mesh/V5_02_GEO.bot","name":"V5_02_GEO","ext":".bot"}
    """
    files = []
    for r, d, f in os.walk(swanFolder):
      for file in f:
        folder=os.path.basename(r)
        if(folder=="results"):
          segPath=r.split('/')
          month=int(segPath[-2])
          year=int(segPath[-3])
          name,ext=os.path.splitext(file)
          files.append({"year":year,"month":month,"path":os.path.join(r, file),"name":name,"ext":ext})
        else:
          name,ext=os.path.splitext(file)
          files.append({"name":name,"ext":ext,"path":os.path.join(r, file)})
    return files
  
  @staticmethod
  def load(filepath,*args,**kwargs):
    """ Loading file
    """
    ext=os.path.splitext(filepath)[1]
    if ext==".bot":return NetCDFSWAN.loadBot(filepath)
    elif ext==".ele":return NetCDFSWAN.loadEle(filepath)
    elif ext==".node":return NetCDFSWAN.loadNode(filepath)
    elif ext==".mat":return NetCDFSWAN.loadMat(filepath)
    elif ext==".json":return NetCDFSWAN.loadJSON(filepath)
    elif ext==".spc":return NetCDFSWAN.loadSpc(filepath,*args,**kwargs)
    else:raise Exception("NetCDFSWAN.load: File type {} is not implemented".format(ext))

  @staticmethod
  def loadJSON(filepath):
    """ Loading json file
    """
    with open(filepath, 'r') as f:
      return json.load(f)
  
  @staticmethod
  def loadBot(filepath):
    """ Loading .bot file
    """
    with open(filepath, 'r') as f:
      array = np.loadtxt(f)
      return array  
  
  @staticmethod
  def loadNode(filepath):
    """ Loading .node file
        Remove first row (only has nnodes) and selecting 1,2 columns(x,y)
    """
    with open(filepath, 'r') as f:
      array = np.loadtxt(f, skiprows=1)
      array = array[:,[1,2]] 
      return array
  
  @staticmethod
  def loadEle(filepath):
    """ Loading .ele file.Remove first row (only has nelem). 
        Ignoring first column (counter) and converting the index instead of id by substracting 1.
    """
    with open(filepath, 'r') as f:
      array = np.loadtxt(f, skiprows=1)
      array = array[:, 1:].astype(int) - 1
      return array
      
  @staticmethod
  def loadMat(filepath,return_datetime=False):
    """ 
    Loading matlab file. 
    Need to convert key to np.datetime64
    Need to sort the keys (datetime) since they might not be in ascending order.
    Need to determine the number of unique datetime and unique variable since the keys might contain multiple variable (e.g. Wind).
    We assume the same number of nodes for each variable.
    
    Returns
    -------
    dict:
      "datetime":np.array((ntime),np.datetime64[s])
      {name of variable}:np.array((ntime,nnode))
      e.g.
      {
        "datetime":np.array(),dtype="datetime64"
        "Windv_x":np.array([[0.,...,178000],[0.,...,178000]]),dtype="f8"
        "Windv_y":np.array([[0.,...,178000],[0.,...,178000]]),dtype="f8"
      }
      
    """
    matfile = loadmat(filepath)
    
    matfile.pop('__header__')
    matfile.pop('__version__')
    matfile.pop('__globals__')
    
    # Extract info from keys only
    timesteps=[]
    for i,key in enumerate(matfile):
      keys=key.split("_")
      name="_".join(keys[:-2]) # To get Windv_x and Windv_y
      date=keys[-2]
      time=keys[-1]
      dt = np.array(datetime(
        int(date[:4]), int(date[4:6]), int(date[6:8]), # Date
        int(time[:2]), int(time[2:4]), int(time[4:6])) # Time
        ,dtype="datetime64[s]")
      timesteps.append({"key":key,"dt":dt,"name":name,"dtkey":date+"_"+time})
    
    # Sort by datetime
    timesteps=sorted(timesteps,key=lambda x:x['dt'])

    # Get shape of output array (ntime,nnode)
    uniqueDT=np.unique(list(map(lambda x:x['dtkey'],timesteps)))
    uniqueNames=np.unique(list(map(lambda x:x['name'],timesteps)))
    shape = (len(uniqueDT),*np.squeeze(matfile[next(iter(matfile))]).shape)
    
    # Initialize dict/array for each variable
    output={}
    for name in uniqueNames:
      output[name]=np.zeros(shape)
    
    # Save array
    for timestep in timesteps:
      name,key,dtkey=timestep['name'],timestep['key'],timestep['dtkey']
      i=np.where(uniqueDT==dtkey)[0][0]
      output[name][i]=np.squeeze(matfile[key])
    
    # Extract datetime from dict
    output['datetime']=np.array(list(map(lambda x:x['dt'],timesteps)))
    return output

  @staticmethod
  def loadSpc(filepath,return_metadata=False,monthOnly=None):
    """
    Aggregates all spc files (stored as .spc instead of .mat) for one month.
    Each datum in the block is multiplied with that block's FACTOR.

    number of timesteps:
      - equal to month, same as mat files.

    number of FACTOR blocks for that timestep:
      - equal to number of lat+lon nodes (e.g. 22 pairs)

    create table for that FACTOR block:
      - columns are each direction "dir" (e.g. 36 across)
      - rows are each frequency "afreq" (e.g. 34 down)
                    
    Parameters
    ----------
    filepath:str
      Path of file
    return_meta:bool
      return metadata only
    
    Returns
    -------
    metadata:
      {
       nlatlng:int
       nfreq:int
       ndir:int
       latlng:np.array((nnode,2))
       freq:np.array((nfreq))
       dir:np.array((ndir))
      }
    output:
    e.g.
      {
        "datetime":np.array(ntime), dtype="datetime64[s]"
        "spectra":np.array((ntime,nnode,nfreq,ndir)), dtype="f8"
      }
    
    """
    re_spcdate = re.compile(r"^20[0-9]{6}[.][0-9]{6}$")  
    
    with open(filepath, "r") as s:
      nlonlat = None
      nfreq = None
      ndir = None
      metadata = True
      while metadata:
        line = s.readline()
        token = line.split()[0]
        if token == "LONLAT":
          nlonlat = int(s.readline().split()[0])
          lonlats=[s.readline().split() for _ in range(nlonlat)]
        elif token=="AFREQ":
          nfreq = int(s.readline().split()[0])
          freq=[s.readline().split() for _ in range(nfreq)]
        elif token=="NDIR":
          ndir = int(s.readline().split()[0])
          dir=[s.readline().split() for _ in range(ndir)]
        elif re.match(re_spcdate, token):
          metadata = False
      
      if(nlonlat is None or nfreq is None or ndir is None):
        raise Exception("Header has a different format")
      
      lonlats=np.array(lonlats).astype("f8")
      freq=np.squeeze(np.array(freq)).astype("f8")
      dir=np.squeeze(np.array(dir)).astype("f8")
      
      if return_metadata:return {"nlatlng":nlonlat,"nfreq":nfreq,"ndir":ndir,"latlng":lonlats,"freq":freq,"dir":dir}
      
      data = True
      output=[]
      datetimes=[]
      
      while data:
        token = line.split()[0]
        if re.match(re_spcdate, token):
          year=int(token[:4]);month=int(token[4:6]);day=int(token[6:8]);hour=int(token[9:11])
          date_time = np.array(datetime(year,month,day,hour),dtype="datetime64[s]")
          array=np.zeros((nlonlat,nfreq,ndir),dtype="f8")
          for i in range(nlonlat):
            FACTOR = s.readline().split()[0].strip()
            if FACTOR == "NODATA":  # no data for this lat-long
              continue
            factor = s.readline().split()[0].strip()
            factor = float(factor)
            array[i]=np.array([np.array([float(int(d)*factor) for d in s.readline().split()]) for afreq in range(nfreq)])
          
          if monthOnly and monthOnly==month:
            datetimes.append(date_time)
            output.append(array)
        else:
            raise Exception(f"date mismatch in load_spc(). reading from stopped")
        line = s.readline()  # ready next line
        if not line:
            data = False
      
      return {"datetime":np.array(datetimes),"spectra":np.array(output)}

  def uploadStatic(self,year=2004,month=1):
    """ Upload static files to S3 (.bot,.ele,.bot, datetime,etc...)
    """
    swanFolder=self.swanFolder
    pbar0=self.pbar0
    pbar=self.pbar
    
    
    if swanFolder is None:raise Exception("swanFolder was not specified")
    _files = NetCDFSWAN.getFiles(swanFolder)
    spcFiles=list(filter(lambda file:"year" in file and "month" in file and file["year"]==year and file["month"]==month and file['ext']=='.spc',_files)) # Get one output folder
    
    if pbar0: pbar0.reset(total=5+1+1+2)
    
    def update(name):
      if pbar0: pbar0.set_description(name);pbar0.update(1)
    
    def update1(name):
      if pbar: pbar.set_description(name);pbar.update(1)

    files=list(filter(lambda file:not "year"in file,_files))
    _bot = list(filter(lambda file:file['ext']==".bot",files))
    _node = list(filter(lambda file:file['ext']==".node",files))
    _ele = list(filter(lambda file:file['ext']==".ele",files))
     
    if len(_bot)>0:self['nodes','bed']=NetCDFSWAN.load(_bot[0]['path']);update("nodes")
    if len(_ele)>0:self['elem','elem']=NetCDFSWAN.load(_ele[0]['path']);update("elem")
    
    if len(_node)>0:
      xy=NetCDFSWAN.load(_node[0]['path'])  
      self['nodes','lon']=xy[:,0];update("lon")
      self['nodes','lat']=xy[:,1];update("lat")
    
    # Upldoad datetime
    self['time','time']=self.datetime;update("time")
    
    if pbar: pbar.reset(total=len(spcFiles))
    # Upldoad stations
    for spcFile in spcFiles:
      name=spcFile['name']
      iIndex=self.stations[name]['start']
      eIndex=self.stations[name]['end']
      id=self.stations[name]['id']
      
      spc=NetCDFSWAN.load(spcFile['path'],return_metadata=True)
      
      latlng=spc['latlng']
      self['snodes','slon',iIndex:eIndex]=latlng[:,0]
      self['snodes','slat',iIndex:eIndex]=latlng[:,1]
      self['snodes','stationid',iIndex:eIndex]=np.zeros(eIndex-iIndex,dtype="int32")+id
      update1(name)
    update('spc')
  
    if pbar: pbar.reset(total=len(self.stations.keys()))
    for name in self.stations:
      id=self.stations[name]['id']
      self['stations','name',id]=name
      update1(name)
    update('stations')
    # Add freq and dir
    self['freq','freq']=spc['freq'];update("freq")
    self['dir','dir']=spc['dir'];update("dir")

    

  def getDatetimeIndex(self,dt):
    """ Find first and last datetime index.
    """
    _datetime=self.datetime
    startDate=dt[0]
    endDate=dt[len(dt)-1]
    startIndex=np.where(_datetime==startDate)[0][0]
    endIndex=np.where(_datetime==endDate)[0][0]+1
    return startIndex,endIndex
    
  def loadRemainingFilestoUpload(self,groupName):
    """ 
    Create temporay json file of all files/groups(called uploadFiles) that needs to be uploaded
    For "s" and "spc" group, "uploadFiles" are the same as the original .mat or .spc files.
    However, for the "t"group, files dont exist explicitly. 
    These "files"/groups needs to be created using the original output files "s". 
    The "uploadFiles" will be grouped using gnode=50
    """
    
    # Get appropriate files (.mat or .spc) for each group ("s","t","spc", "pt")
    if groupName in ["s", "t", "pt"]:
      files=self.matFiles
    elif groupName=="spc":
      files=self.spcFiles
    else: 
      raise Exception("Group-{} not set up".format(groupName))
    
    # Check if json file already exist. 
    # If so, return the json that contains all the files that need to be uploaded.
    path=os.path.join(self.cacheLocation,self.name,groupName+".json")
    if os.path.exists(path):
      uploadFiles=NetCDFSWAN.load(path)
      return files,uploadFiles
    
    # If json file does not exist, create the group/uploadFiles
    # Initialize groups
    if groupName=="spc":
      uploadFiles=files
    elif groupName=="s":
      uploadFiles=list(filter(lambda file: file['name'] not in self.ptList, files))
    elif groupName=="t":
      uploadFiles=list(filter(lambda name: name!="spectra" and name not in self.ptList, list(self.variables.keys())))
      files=list(filter(lambda file: file['path'] not in self.errorlist, files)) # To remove files with errors
    elif groupName=="pt":
      uploadFiles=list(filter(lambda name: name in self.ptList, list(self.variables.keys())))
      files=list(filter(lambda file: file['path'] not in self.errorlist, files)) # To remove files with errors
    else:
      raise Exception("Needs to be s,t,spc,pt")
    
    with open(path,"w+") as f:json.dump(uploadFiles, f)
      
    return files,uploadFiles
  
  def removeUploadedFile(self,groupName,groups):
    """ Remove file from json since it's been uploaded
    """
    if not isinstance(groups,list):groups=[groups] # Make sure it's a list
    path=os.path.join(self.cacheLocation,self.name,groupName+".json")
    with open(path,"w") as f:json.dump(groups, f)
  
  def addErrorFile(self,file):
    """ Add file to the error json file.
    """
    self.errorlist.append(file['path'])
    path = os.path.join(self.folder, "error_list.json")
    with open(path,"w+") as f: f.write(self.errorlist)
  
  def uploadS(self):
    self._uploadSpatial("s")
  def uploadSpc(self):
    self._uploadSpatial("spc")
  def uploadT(self):
    self._uploadTemporal("t")
  def uploadPt(self):
    self._uploadPartition("pt")
  
    
  def _uploadSpatial(self,groupName):
    """ Upload "spatial" output. 
        This is relatively simple since it's uploading directly from Matlab output files to S3.
        Uploading one file at a time.
    """
    showProgress=self.showProgress
    mtname=self.mtname
    variables=self.variables
    stations=self.stations

    _,groups=self.loadRemainingFilestoUpload(groupName)
    
    pbar0=self.pbar0
    if pbar0 is not None: pbar0.reset(total=len(groups))
    
    while groups:
      file=groups.pop(0)
      if pbar0: pbar0.set_description(file['name'])
      if self.logger:self.logger.info("Uploading {}".format(file['path']))      

      try:
        _sub=NetCDFSWAN.load(file['path'],monthOnly=file['month']) # Load matlab or spc file
      except Exception as e:        
        if self.logger: self.logger.info("Couldn't upload {}. Check error_list.json".format(file['path']))
        self.addErrorFile(file)
        if pbar0:pbar0.update(1)
        self.removeUploadedFile(groupName,groups)              
        continue      
        
      dt=_sub.pop("datetime") # Remove datetime from dict since we don't want to upload this
      sIndex,eIndex=self.getDatetimeIndex(dt) # Get datetime index
      
      for key in _sub: # Loop for each variable(e.g Windv_x and Windv_y)
        array=_sub[key]
        name=mtname[key] # Convert matlab/spc name variable to nca name variable
        if groupName=="s":
          self[groupName,name,sIndex:eIndex]=array #upload
        else:
          # Spectral file
          array=np.einsum("abcd->bacd",array) # Need to swap first and second axes
          _sIndex=stations[file['name']]['start']
          _eIndex=stations[file['name']]['end']
          self[groupName,name,_sIndex:_eIndex,sIndex:eIndex]=array # Upload
      
      if pbar0:pbar0.update(1)
      self.removeUploadedFile(groupName,groups)

  
  def _uploadTemporal(self,groupName):
    """ Upload Temporal results
        Gets remaining uploadFiles to upload
        Get node ids for the "group"
        Create memory array to store matlab results for each variable (need to temporary save the results to memory)
        Upload for each variable
    """
    showProgress=self.showProgress
    # mtname=self.mtname
    nnode=self.nnode
    gnode=self.gnode
    ntime=self.ntime
    variables=self.variables

    files,groups=self.loadRemainingFilestoUpload(groupName)

    pbar=self.pbar
    pbar0=self.pbar0
    if pbar0 is not None: pbar0.reset(total=len(groups))
    
    while groups:
      vname=groups.pop(0)
      if pbar0: pbar0.set_description(vname)
      if self.logger is not None:self.logger.info("Uploading group {}".format(vname))
      
      filename=os.path.join(self.cacheLocation,vname+".dat")
      
      variable=variables[vname]
      fileKey=variable["fileKey"]
      _files=list(filter(lambda file:file['name']==fileKey,files))
      
      # Initialize array
      fp = np.memmap(filename, dtype='float32', mode='w+', shape=(nnode,ntime))
      
      if pbar is not None: pbar.reset(total=len(_files))
      # Save matlab info to array
      for _,file in enumerate(_files):
        _sub=NetCDFSWAN.load(file['path'])
        dt=_sub.pop("datetime")
        sIndex,eIndex=self.getDatetimeIndex(dt)
        for key in _sub:
          if key==variable['matfile name']:
            array=_sub[key]
            fp[:,sIndex:eIndex]=array.T
        if pbar:pbar.update(1)
      
      # Save array in memory to S3
      for i in np.arange(0,nnode,gnode):
        _slice=slice(i,np.minimum(nnode,i+gnode))
        self[groupName,vname,_slice]=fp[_slice]
      
      if pbar0:pbar0.update(1)
      self.removeUploadedFile(groupName,groups)

      os.remove(filename)
      



  def _uploadPartition(self,groupName):
    """ Upload Partition results
        Gets remaining uploadFiles to upload
        Get node ids for the "group"
        Create memory array to store matlab results for each variable (need to temporary save the results to memory)
          - (partitions will need additional dimension in the numpy memmap shape)
        Upload for each variable
    """
    showProgress=self.showProgress
    nnode=self.nnode
    gnode=self.gnode
    ntime=self.ntime
    npart=self.npart 
    variables=self.variables

    files,groups=self.loadRemainingFilestoUpload(groupName)
    
    pbar=self.pbar
    pbar0=self.pbar0
    if pbar0 is not None: pbar0.reset(total=len(groups))

    while groups: # for each variable in ["hspt", "drpt", "dspt", "tppt", "stpt", "wlpt"] 
      vname=groups.pop(0)
      if pbar0: pbar0.set_description(vname)
      if self.logger is not None:self.logger.info("Uploading group {}".format(vname))

      filenames = [os.path.join(self.cacheLocation,vname+f"{p:02d}"+".dat") for p in range(1, npart+1)] # create list of file names

      variable=variables[vname]
      fileKey=variable["fileKey"]

      # enter p loop, only loading data into arrays. for each partition file for that variable
      for p in range(npart): 	
	# reset all files needed
	_files=list(filter(lambda file_: (file_['name']==fileKey),files))
        
	# get file name for partition p
	filename=filenames[p]
        
	# grab the memmap file (or create new one if not available)
	fp = np.memmap(filename, dtype='float32', mode='w+', shape=(nnode,ntime))
        
	if pbar is not None: pbar.reset(total=len(_files))
				
        # Save matlab info to array
        for _, file_ in enumerate(_files):
          _sub = NetCDFSWAN.load(file_['path'])
          dt = _sub.pop("datetime")
          sIndex,eIndex=self.getDatetimeIndex(dt)
          for key in _sub:                               # because more than one partition in the file "HsPT01_", "HsPT02_", etc.
            key_name, key_num = key[:-2], int(key[-2:])  # split "HsPT01" into "HsPT" and 1
            if key_num == p+1:                           # check key name. no key "0", starts from "1". 
              array = _sub[key]                          # e.g. "HsPT01" data array
              fp[:, sIndex:eIndex]=array.T               # load array into memmap
          if pbar:pbar.update(1)
					
      if pbar is not None: pbar.reset(total=nnode/gnode) # e.g. total is 595 if gnode == 300
			
      # enter p loop again, this time for uploading
      for p in range(npart):     
        filename = filenames[p] # get file name for partition p
        self._uploadPartitionSliceFile(groupName, vname, filename, p) # enter g loop (see below)
        if pbar:pbar.update(1)
        os.remove(filename)
      
      if pbar0:pbar0.update(1)
      self.removeUploadedFile(groupName,groups)

  
  def _uploadPartitionSliceFile(self, groupName, vname, filename, p):
    nnode=self.nnode
    gnode=self.gnode
    ntime=self.ntime
    pbar=self.pbar
    
    # read the memmap file for partition p
    fp = np.memmap(filename, dtype='float32', mode='r', shape=(nnode,ntime))
    
    # g loop (Save array in memory to S3)
    for i in np.arange(0,nnode,gnode):
      _slice=slice(i,np.minimum(nnode,i+gnode))
      self[groupName, vname, _slice, p] = fp[_slice]
    
