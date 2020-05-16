import os
import json
from netcdfswan import NetCDFSWAN
import numpy as np
import logging

from dataTest import elem,time,lat,lon,bed,slat,slon,freq,dir,spcgroup,variables

def test_NetCDFSWAN_write():
  swanFolder='./output'
  jsonFile='./json/demo.json'
  input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder,year=2000,month=1)
  
  swan=NetCDFSWAN(input)
  # Write
  swan.uploadStatic(year=2000)
  swan.uploadS()
  swan.uploadT()
  swan.uploadSpc()

def test_NetCDFSWAN():

  
  input={
    "name":"test1",
    "bucket":"uvic-bcwave",
    "cacheLocation":"../../s3",
    "localOnly":False
  }
  
  swan=NetCDFSWAN(input)

  # Read
  np.testing.assert_array_equal(swan["nodes","bed"], bed)
  np.testing.assert_array_equal(swan["elem","elem"], elem)
  np.testing.assert_array_equal(swan["time","time"], time)
  np.testing.assert_array_equal(swan["nodes","lat"], lat)
  np.testing.assert_array_equal(swan["nodes","lon"], lon)
  
  
  
  np.testing.assert_array_equal(swan["freq","freq"], freq)
  np.testing.assert_array_equal(swan["dir","dir"], dir)

  np.testing.assert_array_equal(swan["s","u10"], variables['WIND']['Windv_x'])
  np.testing.assert_array_equal(swan["s","v10"], variables['WIND']['Windv_y'])
  np.testing.assert_array_equal(swan["s","hs"], variables['HS']['Hsig'])
  
  # TODO: Replace name of variables below
  # np.testing.assert_array_equal(swan["s","tps"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","tmm10"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","tm01"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","tm02"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","pdir"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","dir"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","dspr"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","qp"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","transpx"], variables['WIND']['Windv_x'])
  # # np.testing.assert_array_equal(swan["s","transpy"], variables['WIND']['Windv_x'])
  
  
  np.testing.assert_array_equal(swan["t","u10"], variables['WIND']['Windv_x'].T)
  np.testing.assert_array_equal(swan["t","v10"], variables['WIND']['Windv_y'].T)
  np.testing.assert_array_equal(swan["t","hs"], variables['HS']['Hsig'].T)
  
  # TODO: Replace name of variables below
  # np.testing.assert_array_equal(swan["t","tps"], variables['tps'].T)
  # np.testing.assert_array_equal(swan["t","tmm10"], variables['tmm10'].T)
  # np.testing.assert_array_equal(swan["t","tm01"], variables['tm01'].T)
  # np.testing.assert_array_equal(swan["t","tm02"], variables['tm02'].T)
  # np.testing.assert_array_equal(swan["t","pdir"], variables['pdir'].T)
  # np.testing.assert_array_equal(swan["t","dir"], variables['dir'].T)
  # np.testing.assert_array_equal(swan["t","dspr"], variables['dspr'].T)
  # np.testing.assert_array_equal(swan["t","qp"], variables['qp'].T)
  # np.testing.assert_array_equal(swan["t","transpx"], variables['transpx'].T)
  # np.testing.assert_array_equal(swan["t","transpy"], variables['transpy'].T)  
  
  
  # TODO: need to change default values in spcgroup
  # TODO: swan["spc","spectra"] does not work...I'll have to check s3-netcdf
  
  # print(swan.groups['spc'].child)
  # for i,id in enumerate(stations):
  #   nsnode =stations[id]
  #   for j in range(nsnode):
  #     np.testing.assert_array_equal(swan["spc","spectra",i,j], spcgroup['spectra'][0,1])
  # np.testing.assert_array_equal(swan["spc","spectra",0,1], spcgroup['spectra'][0,1])
  
  # x=swan["spc","spectra",1,1]
  
  # np.testing.assert_array_equal(swan["spc","spectra",1,0], spcgroup['spectra'][1,0])
  # np.testing.assert_array_equal(swan["spc","spectra",8], spcgroup['spectra'][8])



def test_NetCDFSWAN_logger():
  logging.basicConfig(
        filename=os.path.join('./data',"progress.log"),
        level=logging.DEBUG,
        format="%(levelname)s %(asctime)s  --| %(message)s"
        )
  logger = logging.getLogger()
  
  try:
    swanFolder='./data'
    jsonFile='./json/demo.json'
    input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder,year=2000,month=1)
    swan=NetCDFSWAN(input,logger=logger)
    swan.uploadStatic()
    swan.uploadS()
    swan.uploadT()
    swan.uploadSpc()
  except Exception as err:
    logger.error(err)
  
  
if __name__ == "__main__":
  # test_NetCDFSWAN_write()
  test_NetCDFSWAN()
  