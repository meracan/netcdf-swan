from datetime import datetime
from netCDF4 import num2date, date2num
from netcdfswan import NetCDFSWAN
import random
from scipy.io import loadmat

def test_uvicswan():
  
  # this will need the location of the swan data on the server 
  # (relative to test_netcdfswan, or use absolute path)
  data_location = ?
  
  # create input json object (needs nca object)
  # make sure localOnly is set to False
  inputFile = "./BCSWANv5.json"
  with open(inputFile, "r") as f:
    inputJson = json.load(f)
  
  # load from s3 bucket
  swan = NetCDFSWAN(inputJson)
  
  # constants
  u, c = "hours since 1970-01-01 00:00:00.0", "gregorian"
  tSTART = int(date2num(datetime(2004,1,1,0,0), units=u, calendar=c))
  
  # trial list.
  mats = {
    "u10":("WIND", "Windv_x"),
    "v10":("WIND", "Windv_y"),
    "tps":("TPS", "TPsmoo"),
    "tm01":("TM01", "Tm01"),
    "dir":("DIR", "Dir")
  }
  for mat in mats.items():
    var, mVAR, mvar = mat[0], mat[1][0], mat[1][1]
    
    # create 5 random timesteps to check
    for i in range(5):
      y = random.randrange(2004, 2017)
      m = random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      
      mfilepath = data_location + str(y)+"/"+f"{m:02d}"+"/results/"+mVAR+".mat"
      try: matdata = loadmat(mfilepath)
      except: print(f"couldnt read {mfilepath}"); continue
      
      key = mvar+"_"+str(y)+f"{m:02d}"+f"{d:02d}"+"_"+f"{h:02d}"+"0000"
      local_nodes = matdata[key][0]
      rmote_nodes = swan["s", var, t]
      
      np.testing.assert_array_equal(local_nodes, rmote_nodes)
      print(f"{key}: t={t} OK")