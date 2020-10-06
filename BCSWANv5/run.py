import os
import json
from netcdfswan import NetCDFSWAN


if __name__ == "__main__":
  # swanFolder='../s3/data'
  # jsonFile='BCSWANv5/BCSWANv5.json'
  # input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder,year=2004,month=1)
  
  
  
  input={
    "name":"SWANv5a",
    "swanFolder":'../s3/data',
    "bucket":"uvic-bcwave",
    "showProgress":True,
    "cacheLocation":"../s3",
    "localOnly":False
  }
  
  
  swan=NetCDFSWAN(input)
  # print(swan.groups['t'].child)

  # swan.s3.delete()
  # print(swan.groups['s'].child)
  # Write
  # swan.uploadStatic(year=2004)
  swan.uploadS()
  # swan.uploadT()
  # swan.uploadSpc()  