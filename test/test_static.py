import os
import json
from netcdfswan import NetCDFSWAN

tmpFolder="../s3/swandata"


def test_getFiles():
  print(NetCDFSWAN.getFiles(tmpFolder))

def test_load():
  print(NetCDFSWAN.load(os.path.join(tmpFolder,"2000/1/results/HS.mat"))['Hsig'].shape)  
  print(NetCDFSWAN.load(os.path.join(tmpFolder,"2000/1/results/WIND.mat")))  
  print(NetCDFSWAN.load(os.path.join(tmpFolder,"Mesh/dummy.bot")).shape)
  print(NetCDFSWAN.load(os.path.join(tmpFolder,"Mesh/dummy.ele")).shape)
  
def test_printInfo():
  NetCDFSWAN.printMatKeys(tmpFolder,year=2000)
  NetCDFSWAN.printSpcShape(tmpFolder,year=2000)

def test_getSpectralStationMetadata():
  obj=NetCDFSWAN.getSpectralStationMetadata(tmpFolder,year=2000)
  with open(os.path.join(tmpFolder,"demo.spectral.json"),"w+") as f:
    json.dump(obj,f,indent=2)

def test_preparingFunctions():
  obj=NetCDFSWAN.prepareInputJSON('./test/json/demo.json',tmpFolder,year=2000)
  with open(os.path.join(tmpFolder,"demo.prepared.json"),"w+") as f:
    json.dump(obj,f,indent=2)


if __name__ == "__main__":
  # test_getFiles()
  # test_load()
  # test_printInfo()
  # test_getSpectralStationMetadata()
  test_preparingFunctions()