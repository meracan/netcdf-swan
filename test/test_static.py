import os
import json
from netcdfswan import NetCDFSWAN

def test_getFiles():
  print(NetCDFSWAN.getFiles("./data"))

def test_load():
  print(NetCDFSWAN.load("./data/2004/01/results/HS_6hr.mat")['Hsig'].shape)  
  print(NetCDFSWAN.load("./data/2004/01/results/WIND_6hr.mat"))  
  print(NetCDFSWAN.load("./data/Mesh/V5_02_GEO.bot").shape)
  print(NetCDFSWAN.load("./data/Mesh/V5_02_GEO.ele").shape)
  
def test_printInfo():
  NetCDFSWAN.printMatKeys("./data")
  NetCDFSWAN.printSpcShape("./data")

def test_getSpectralStationMetadata():
  obj=NetCDFSWAN.getSpectralStationMetadata("./data")
  with open("./output/demo.spectral.json","w+") as f:
    json.dump(obj,f,indent=2)

def test_preparingFunctions():
  obj=NetCDFSWAN.prepareInputJSON('./json/demo.json',"./data")
  with open("./output/demo.prepared.json","w+") as f:
    json.dump(obj,f,indent=2)



if __name__ == "__main__":
  test_getFiles()
  test_load()
  test_printInfo()
  test_getSpectralStationMetadata()
  test_preparingFunctions()