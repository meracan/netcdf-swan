import os
import numpy as np
# import scipy
from scipy.io import savemat
from dataTest import elem,time,lat,lon,bed,slat,slon,freq,dir,variables,spcgroup,stations,nodes
from netcdfswan import NetCDFSWAN
from datetime import datetime,timedelta
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

  #pp.pprint(mdict)
  #print(os.path.abspath(filePath))
  savemat(filePath, mdict)


  
  
def create_spc(filePath,dic,station):
  # dic:{"datetime":datetime64,"datetime":np.array((ntime,nsnode,nfreq,ndir))}
  freq=dic['freq']
  dir=dic['dir']
  dt=dic['datetime']
  spectra=dic['spectra']
  latlng=station['latlng']
  stationId=station['id']
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
    
    for i,_ in enumerate(dt):
      dtStr=_.astype(object).strftime("%Y%m%d.%H%M%S")
      f.write("{}                         date and time\n".format(dtStr))
      for inode,_ in enumerate(latlng):
        # print(stationId,inode,i)
        f.write("FACTOR\n")
        factor=1
        array=(spectra[stationId,inode,i]/factor).astype("i4")
        
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
  folder="./output"
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
    print(date)
    if day==1 and hour ==0:
     
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
        nsnodes=station['nsnodes']
        dic={
            "datetime":dt,
            "freq":freq,
            "dir":dir,
            "spectra":spcgroup['spectra'][:,:,startIndex:endIndex]
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
  # TODO: read matlab file,spc and make sure it's the same array
  # np.testing.assert_array_equal(NetCDFSWAN.load('./output/2000/12/results/HS.mat',variables['HS'][['Hsig']])
  None
        
if __name__ == "__main__":
  create_data()
  # check_data()
  