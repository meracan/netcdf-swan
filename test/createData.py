import os
import numpy as np
# import scipy
from scipy.io import savemat
from dataTest import elem,time,lat,lon,bed,slat,slon,freq,dir,variables,spcgroup,stations,nodes
from netcdfswan import NetCDFSWAN
from datetime import datetime,timedelta


tmpFolder="../s3/swandata"

def create_bot(filePath,array):
  np.savetxt(filePath, array, delimiter=' ')   

def create_ele(filePath,array):
  index=np.arange(1,len(array)+1)
  array=array+1
  array=np.column_stack((index,array))
  np.savetxt(filePath, array,header='{} 3 0'.format(len(array)))

def create_node(filePath,array):
  array=np.column_stack((np.arange(1,len(array)+1),array,array,np.ones(len(array))))
  np.savetxt(filePath, array,header='{} 2 0 1'.format(len(array)))

def create_mat(filePath,dic):
  # dic:{"datetime":np.array(ntime,datetime64),"hs":np.array((ntime,nnode))}
  dts=dic.pop('datetime')
  mdict={}
  for i,dt in enumerate(dts):
    dtStr=dt.astype(object).strftime("%Y%m%d_%H%M%S")
    for name in dic:
      key="{}_{}".format(name,dtStr)
      mdict[key]=dic[name][i]
      
  mdict['__header__']="dummy header"
  mdict['__version__']="dummy version"
  mdict['__globals__']="dummy globals"


  savemat(filePath, mdict)


  
  
def create_spc(filePath,dic,station):
  # dic:{"datetime":datetime64,"datetime":np.array((ntime,nsnode,nfreq,ndir))}
  freq=dic['freq']
  dir=dic['dir']
  dt=dic['datetime']
  spectra=dic['spectra']
  latlng=station['latlng']

  with open(filePath,"w+") as f:
    f.write("SWAN   1                                Swan standard spectral file, version\n")
    f.write("$   Data produced by SWAN version 41.31    \n")
    f.write("$   Project: WCWI_V5         ;  run number: 01  \n")
    f.write("TIME                                    time-dependent data\n")
    f.write("     1                                  time coding option\n")
    f.write("LONLAT                                  locations in spherical coordinates\n")
    f.write("     {}                                  number of locations\n".format(len(latlng)))
    
    arrayStr = np.array2string(latlng,separator=' ').replace('[',"").replace(']',"")
    f.write("{}\n".format(arrayStr))

    f.write("AFREQ absolute frequencies in Hz\n")
    f.write("{} number of frequencies\n".format(len(freq)))
    arrayStr = np.array2string(freq,separator='\n').replace('[',"").replace(']',"")
    f.write("{}\n".format(arrayStr))

    f.write("NDIR                                    spectral nautical directions in degr\n")
    f.write("   {}                                    number of directions\n".format(len(dir)))
    arrayStr = np.array2string(dir,separator='\n').replace('[',"").replace(']',"")
    f.write("{}\n".format(arrayStr))
    
    f.write("QUANT\n")
    f.write("    1                                   number of quantities in table\n")
    f.write("VaDens                                  variance densities in m2/Hz/degr\n")
    f.write("m2/Hz/degr                              unit\n")
    f.write("-0.9900E+02                             exception value\n")
    
    # ------- To test spinup, add extra date --------------
    if dt[0].astype('datetime64[M]').astype(int) % 12 + 1==1:
      _d = dt[0]-np.timedelta64(10,"D")
      
      dtStr=_d.astype(object).strftime("%Y%m%d.%H%M%S")
      f.write("{}                         date and time\n".format(dtStr))
      for inode,_ in enumerate(latlng):
        # print(stationId,inode,i)
        f.write("FACTOR\n")
        factor=1
        array=(spectra[inode,0]/factor).astype("i4")
        
        arrayStr = np.array2string(array,separator=',').replace(" ","").replace('[',"").replace(']',"").replace(","," ")
        f.write("{}\n".format(factor))
        f.write("{}\n".format(arrayStr))
    # -------------------------------------------------------
    for i,_ in enumerate(dt):
      dtStr=_.astype(object).strftime("%Y%m%d.%H%M%S")
      f.write("{}                         date and time\n".format(dtStr))
      for inode,_ in enumerate(latlng):
        # print(stationId,inode,i)
        f.write("FACTOR\n")
        factor=1
        array=(spectra[inode,i]/factor).astype("i4")
        
        arrayStr = np.array2string(array,separator=',').replace(" ","").replace('[',"").replace(']',"").replace(","," ")
        f.write("{}\n".format(factor))
        f.write("{}\n".format(arrayStr))
  
def getDatetimeIndex(_all,dt):
  """ Find first and last datetime index.
  """
  startDate=dt[0]
  endDate=dt[len(dt)-1]
  
  startIndex=np.where(_all==startDate)[0][0]
  if(endDate>_all[len(_all)-1]):
     endIndex=len(_all)+1
  else:
    endIndex=np.where(_all==endDate)[0][0]+1
   
  return startIndex,endIndex
  
  
def create_folders(folder,year,month):
  # print(folder,year)
  yearFolder=os.path.join(folder,str(year))
  if not os.path.exists(yearFolder):os.mkdir(yearFolder)
  monthFolder=os.path.join(yearFolder,str(month))
  if not os.path.exists(monthFolder):os.mkdir(monthFolder)
  resultsFolder=os.path.join(monthFolder,"results")
  if not os.path.exists(resultsFolder):os.mkdir(resultsFolder)
  return resultsFolder

def create_data():
  folder=tmpFolder
  if not os.path.exists(folder):os.mkdir(folder)
  meshFolder=os.path.join(folder,"Mesh")
  if not os.path.exists(meshFolder):os.mkdir(meshFolder)
  
  # Create .bot,.ele,.bot
  create_ele(os.path.join(meshFolder,"dummy.ele"),elem)
  create_bot(os.path.join(meshFolder,"dummy.bot"),bed)
  create_node(os.path.join(meshFolder,"dummy.node"),nodes)
  
  
  # Create output
  for iday,t in enumerate(time):
    date=t.astype(object)
    year = date.year
    month = date.month
    day=date.day
    hour=date.hour
    if day==1 and hour ==0:
      print(date, "...")
      if month+1>12:
        endDate= datetime(year+1,1,1)
      else:
        endDate= datetime(year,month+1,1)
      
      if np.datetime64(endDate)>time[len(time)-1]:
         endDate=(time[len(time)-1]).astype(object)
    
      dt=np.arange(date,endDate+timedelta(hours=1), timedelta(hours=1)).astype("datetime64[s]")
      
      startIndex,endIndex=getDatetimeIndex(time,dt)
      resultsFolder=create_folders(folder,year,month)
      
      for i,s in enumerate(stations):
        station=stations[s]
        sIndex=station['start']
        eIndex=station['end']
        nsnodes=station['nsnodes']
        dic={
            "datetime":dt,
            "freq":freq,
            "dir":dir,
            "spectra":spcgroup['spectra'][sIndex:eIndex,startIndex:endIndex]
        }
        
        create_spc(os.path.join(resultsFolder,"{}.spc".format(s)),dic,station)
      for name in variables:
        variable=variables[name]
        filePath=os.path.join(resultsFolder,name+".mat")
        dic={}
        for _ in variable:
          dic[_]=variable[_][startIndex:endIndex]
        
        dic['datetime']=dt
        create_mat(filePath,dic)
        
def check_data():

  for mkey in variables.keys():
    for mmkey in list(variables[mkey].keys()):
      for month in range(1, 13):
        NCDFS = np.array(NetCDFSWAN.load(os.path.join(tmpFolder,f'2000/{str(month)}/results/{mkey}.mat'))[mmkey])
        start = int(NCDFS[0][0]//10) # actual value is similar to the index
        end = start + NCDFS.shape[0]
        v = variables[mkey][mmkey][start: end]  # 0 - 745
        np.testing.assert_array_equal(NCDFS, v)
      print(f"{mkey}.mat ok")


  for i, station in enumerate(stations):
    n = stations[station]["nsnodes"]
    sts = 0
    for month in range(1, 13):
      NCDFS = NetCDFSWAN.load(os.path.join(tmpFolder,f'2000/{str(month)}/results/{station}.spc'))["spectra"]

      ts = NCDFS.shape[0]
      spg = spcgroup["spectra"][i, :n]

      for node in range(n):
        spg_n = spg[node][sts:sts+ts]
        NCDFS_n = NCDFS[:, node]
        np.testing.assert_array_equal(spg_n, NCDFS_n)
        #print(f"station {station} node {node}:", NCDFS_n.shape, spg_n.shape, " start index:", sts)
      sts += ts - 1
    print(f"{station}.spc ok")

        
if __name__ == "__main__":
  create_data()
  # check_data()