import os
import json
from netcdfswan import NetCDFSWAN

# year=2000 is specified since demo starts in 2000. 

def test_getFiles():
  print(NetCDFSWAN.getFiles("./data"))

def test_load():
  print(NetCDFSWAN.load("./output/2000/1/results/HS.mat")['Hsig'].shape)  
  print(NetCDFSWAN.load("./output/2000/1/results/WIND.mat"))  
  print(NetCDFSWAN.load("./output/Mesh/dummy.bot").shape)
  print(NetCDFSWAN.load("./output/Mesh/dummy.ele").shape)
  
def test_printInfo():
  NetCDFSWAN.printMatKeys("./output",year=2000)
  NetCDFSWAN.printSpcShape("./output",year=2000)

def test_getSpectralStationMetadata():
  obj=NetCDFSWAN.getSpectralStationMetadata("./output",year=2000)
  with open("./output/demo.spectral.json","w+") as f:
    json.dump(obj,f,indent=2)

def test_preparingFunctions():
  obj=NetCDFSWAN.prepareInputJSON('./json/demo.json',"./output",year=2000)
  with open("./output/demo.prepared.json","w+") as f:
    json.dump(obj,f,indent=2)


if __name__ == "__main__":
  test_getFiles()
  test_load()
  test_printInfo()
  test_getSpectralStationMetadata()
  test_preparingFunctions()