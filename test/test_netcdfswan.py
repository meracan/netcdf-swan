import os
import json
from netcdfswan import NetCDFSWAN
import numpy as np
import logging

from .dataTest import elem,time,lat,lon,bed,slat,slon,freq,dir,sgroup,spcgroup


def test_NetCDFSWAN():
  swanFolder='./data'
  jsonFile='./json/demo.json'
  input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder)
  swan=NetCDFSWAN(input)
  # Write
  swan.uploadStatic()
  swan.uploadS()
  swan.uploadT()
  swan.uploadSpc()
  
  # Read
  np.testing.assert_array_equal(swan["elem","elem"], elem)
  np.testing.assert_array_equal(swan["time","time"], time)
  np.testing.assert_array_equal(swan["nodes","lat"], lat)
  np.testing.assert_array_equal(swan["nodes","lon"], lon)
  np.testing.assert_array_equal(swan["nodes","bed"], bed)
  np.testing.assert_array_equal(swan["stations","lat"], slat)
  np.testing.assert_array_equal(swan["stations","long"], slon)
  np.testing.assert_array_equal(swan["freq","freq"], freq)
  np.testing.assert_array_equal(swan["dir","dir"], dir)
  
  np.testing.assert_array_equal(swan["s","u10"], sgroup['u10'])
  np.testing.assert_array_equal(swan["s","v10"], sgroup['v10'])
  np.testing.assert_array_equal(swan["s","hs"], sgroup['hs'])
  np.testing.assert_array_equal(swan["s","tps"], sgroup['tps'])
  np.testing.assert_array_equal(swan["s","tmm10"], sgroup['tmm10'])
  np.testing.assert_array_equal(swan["s","tm01"], sgroup['tm01'])
  np.testing.assert_array_equal(swan["s","tm02"], sgroup['tm02'])
  np.testing.assert_array_equal(swan["s","pdir"], sgroup['pdir'])
  np.testing.assert_array_equal(swan["s","dir"], sgroup['dir'])
  np.testing.assert_array_equal(swan["s","dspr"], sgroup['dspr'])
  np.testing.assert_array_equal(swan["s","qp"], sgroup['qp'])
  np.testing.assert_array_equal(swan["s","transpx"], sgroup['transpx'])
  np.testing.assert_array_equal(swan["s","transpy"], sgroup['transpy'])
  
  np.testing.assert_array_equal(swan["t","u10"], sgroup['u10'].T)
  np.testing.assert_array_equal(swan["t","v10"], sgroup['v10'].T)
  np.testing.assert_array_equal(swan["t","hs"], sgroup['hs'].T)
  np.testing.assert_array_equal(swan["t","tps"], sgroup['tps'].T)
  np.testing.assert_array_equal(swan["t","tmm10"], sgroup['tmm10'].T)
  np.testing.assert_array_equal(swan["t","tm01"], sgroup['tm01'].T)
  np.testing.assert_array_equal(swan["t","tm02"], sgroup['tm02'].T)
  np.testing.assert_array_equal(swan["t","pdir"], sgroup['pdir'].T)
  np.testing.assert_array_equal(swan["t","dir"], sgroup['dir'].T)
  np.testing.assert_array_equal(swan["t","dspr"], sgroup['dspr'].T)
  np.testing.assert_array_equal(swan["t","qp"], sgroup['qp'].T)
  np.testing.assert_array_equal(swan["t","transpx"], sgroup['transpx'].T)
  np.testing.assert_array_equal(swan["t","transpy"], sgroup['transpy'].T)  
  
  np.testing.assert_array_equal(swan["spc","spectra"], spcgroup['spectra'])  

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
    input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder)
    swan=NetCDFSWAN(input,logger=logger)
    
    swan.uploadStatic()
    swan.uploadS()
    swan.uploadT()
    swan.uploadSpc()
  except Exception as err:
    logger.error(err)
  
  
if __name__ == "__main__":
  test_NetCDFSWAN_logger()
  