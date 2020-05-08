import os
import json
from netcdfswan import NetCDFSWAN

swanFolder='../../s3/data'

def test_getFiles():
  print(NetCDFSWAN.getFiles(swanFolder))

def test_printInfo():
  NetCDFSWAN.printMatKeys(swanFolder)
  NetCDFSWAN.printSpcShape(swanFolder)

def test_getSpectralStationMetadata():
  obj=NetCDFSWAN.getSpectralStationMetadata(swanFolder)
  with open(os.path.join(swanFolder,"spectral.json"),"w+") as f:
    json.dump(obj,f,indent=2)

def test_preparingFunctions():
  obj=NetCDFSWAN.prepareInputJSON('BCSWANv5.json',swanFolder)
  with open(os.path.join(swanFolder,"intput.json"),"w+") as f:
    json.dump(obj,f,indent=2)



if __name__ == "__main__":
  # test_getFiles()
  # test_printInfo()
  # test_getSpectralStationMetadata()
  # test_preparingFunctions()