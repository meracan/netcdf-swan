import os
import json
from netcdfswan import NetCDFSWAN


if __name__ == "__main__":
  swanFolder='../../s3/data'
  jsonFile='BCSWANv5.json'
  input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder,year=2004,month=1)
  
  swan=NetCDFSWAN(input)

  # swan.s3.delete()
  # print(swan.groups['s'].child)
  # Write
  # swan.uploadStatic(year=2004)
  swan.uploadS()
  # swan.uploadT()
  # swan.uploadSpc()  