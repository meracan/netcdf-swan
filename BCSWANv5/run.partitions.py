import os
import json
from netcdfswan import NetCDFSWAN
import numpy as np

if __name__ == "__main__":
  swanFolder='../data'
  jsonFile='BCSWANv6/BCSWANv6.json'
  
  
  
  
  
  
  input={
    "name":"SWANv6",
    "swanFolder":'../data',
    "bucket":"uvic-bcwave",
    "showProgress":True,
    "memorySize":40,
    "cacheSize":100,
    "cacheLocation":"../data",
    "localOnly":False
  }
  with NetCDFSWAN(input) as swan:
    swan.uploadPt("hspt")
    # ntime     = swan.obj['dimensions'].get('ntime')
    # startDate = swan.obj['metadata'].get('startDate')
    # timeStep  = swan.obj['metadata'].get('timeStep(h)')
    # startDate=np.datetime64(startDate)
    # datetime  = startDate+np.arange(ntime)*np.timedelta64(timeStep, 'h')
    # swan['time','time']=datetime
    # print(datetime)
    # print(swan['time','time'])
    
    # dirbin=np.array([265,255,245,235,225,215,205,195,185,175,165,155,145,135,125,115,105, 95, 85, 75, 65, 55, 45, 35, 25, 15,  5, -5,-15,-25,-35,-45,-55,-65,-75,-85]) 
    # swan['dirbin','dirbin']=dirbin
    
    
    # swan.uploadStatic(year=2004)
    # swan.uploadS()
    # swan.uploadT()
    # swan.uploadSpc()
  
  