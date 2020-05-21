import os
import json
from netcdfswan import NetCDFSWAN
import numpy as np
import logging

from dataTest import elem,time,lat,lon,bed,slat,slon,freq,dir,spcgroup,variables, stations

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
  np.testing.assert_array_equal(swan["s","tps"], variables['TPS']['TPsmoo'])
  np.testing.assert_array_equal(swan["s","tmm10"], variables['TMM10']['Tm_10'])
  np.testing.assert_array_equal(swan["s","tm01"], variables['TM01']['Tm01'])
  np.testing.assert_array_equal(swan["s","tm02"], variables['TM02']['Tm02'])
  np.testing.assert_array_equal(swan["s","pdir"], variables['PDIR']['Pdir'])
  np.testing.assert_array_equal(swan["s","dir"], variables['DIR']['Dir'])
  np.testing.assert_array_equal(swan["s","dspr"], variables['DSPR']['Dspr'])
  np.testing.assert_array_equal(swan["s","qp"], variables['QP']['Qp'])
  np.testing.assert_array_equal(swan["s","transpx"], variables['TRANSP']['Transp_x'])
  np.testing.assert_array_equal(swan["s","transpy"], variables['TRANSP']['Transp_y'])

  np.testing.assert_array_equal(swan["t","u10"], variables['WIND']['Windv_x'].T)
  np.testing.assert_array_equal(swan["t","v10"], variables['WIND']['Windv_y'].T)
  np.testing.assert_array_equal(swan["t","hs"], variables['HS']['Hsig'].T)

  np.testing.assert_array_equal(swan["t","tps"], variables['TPS']['TPsmoo'].T)
  np.testing.assert_array_equal(swan["t","tmm10"], variables['TMM10']['Tm_10'].T)
  np.testing.assert_array_equal(swan["t","tm01"], variables['TM01']['Tm01'].T)
  np.testing.assert_array_equal(swan["t","tm02"], variables['TM02']['Tm02'].T)
  np.testing.assert_array_equal(swan["t","pdir"], variables['PDIR']['Pdir'].T)
  np.testing.assert_array_equal(swan["t","dir"], variables['DIR']['Dir'].T)
  np.testing.assert_array_equal(swan["t","dspr"], variables['DSPR']['Dspr'].T)
  np.testing.assert_array_equal(swan["t","qp"], variables['QP']['Qp'].T)
  np.testing.assert_array_equal(swan["t","transpx"], variables['TRANSP']['Transp_x'].T)
  np.testing.assert_array_equal(swan["t","transpy"], variables['TRANSP']['Transp_y'].T)

  for i, station in enumerate(stations):
    n = stations[station]["nsnodes"]
    swn = swan["spc", "spectra", i, :n]
    spg = spcgroup["spectra"][i, :n]
    if n > 1:
      for node in range(n):
        swn_n = swn[node]
        spg_n = spg[node]
        np.testing.assert_array_equal(swn_n, spg_n)
        print(f"station {station} node {node} match")
    else:
      spg_n = spg[0]
      np.testing.assert_array_equal(swn, spg_n)
      print(f"station {station} node 0 match")


def test_NetCDFSWAN_logger():
  logging.basicConfig(
        filename=os.path.join('./data',"progress.log"),
        level=logging.DEBUG,
        format="%(levelname)s %(asctime)s   %(message)s"
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
  test_NetCDFSWAN_write()
  test_NetCDFSWAN()
  