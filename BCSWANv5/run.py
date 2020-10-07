import os
import json
from netcdfswan import NetCDFSWAN


if __name__ == "__main__":
  # swanFolder='../s3/data'
  # jsonFile='BCSWANv5/BCSWANv5.json'
  # input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder,year=2004,month=1)
  
  
  
  input={
    "name":"SWANv5",
    "swanFolder":'../s3/data',
    "bucket":"uvic-bcwave",
    "showProgress":True,
    "memorySize":10,
    "cacheSize":100,
    "cacheLocation":"../s3",
    "localOnly":False
  }
  
  
  swan=NetCDFSWAN(input)
  # swan.uploadStatic(year=2004)
  # swan.uploadS()
  # swan.uploadT()
  # swan.uploadSpc()
  
  