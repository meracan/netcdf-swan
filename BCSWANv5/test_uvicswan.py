from datetime import datetime
from netCDF4 import num2date, date2num
from netcdfswan import NetCDFSWAN
import random
import numpy as np
from scipy.io import loadmat
import json

# localOnly should be set to False
# this will need the location of the swan data on the server 
# (relative to test_netcdfswan, or use absolute path)
data_location = ?
  
# constants
u, c = "hours since 1970-01-01 00:00:00.0", "gregorian"
tSTART = int(date2num(datetime(2004,1,1,0,0), units=u, calendar=c))  
num_stations_to_test = 13
num_timesteps_to_test = 3
  
# create input json object. inputFile path depends where this script is run
#inputFile = "./BCSWANv5.json"
inputFile = "./netcdf-swan/BCSWANv5/BCSWANv5.json" # set localOnly to False
with open(inputFile, "r") as f:
  inputJson = json.load(f)
inputJson["showProgress"] = False
  
def test_uvicswan_mat():
  global num_timesteps_to_test
  
  # load from s3 bucket
  swan = NetCDFSWAN(inputJson)
  
  # trial list. Only chooses a few of them
  mats = {
    "u10":("WIND", "Windv_x"),
    "v10":("WIND", "Windv_y"),
    "tps":("TPS", "TPsmoo"),
    "tm01":("TM01", "Tm01"),
    "dir":("DIR", "Dir")
  }
  for mat in mats.items():
    var, mVAR, mvar = mat[0], mat[1][0], mat[1][1]
    
    # create random timesteps to check
    for i in range(num_timesteps_to_test):
      y = random.randrange(2004, 2017)
      m = random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      dateprint = num2date(t+tSTART, units=u, calendar=c)
      
      mfilepath = data_location +"/"+ str(y)+"/"+f"{m:02d}"+"/results/"+mVAR+".mat"
      try: matdata = loadmat(mfilepath)
      except: print(f"couldnt read {mfilepath}"); continue
      
      key = mvar+"_"+str(y)+f"{m:02d}"+f"{d:02d}"+"_"+f"{h:02d}"+"0000"
      local_nodes = matdata[key][0]
      rmote_nodes = swan["s", var, t][0]
      
      np.testing.assert_array_equal(local_nodes, rmote_nodes)
      print(f"{key} {dateprint} OK")


def test_uvicswan_spc():
  global num_stations_to_test, num_timesteps_to_test
  
  # load from s3 bucket
  swan = NetCDFSWAN(inputJson)
  
  stations = json.loads(swan.info()["metadata"]["stations"])
  
  # check stations
  for station_name, values in stations.items():
    if num_stations_to_test <= 0: break
    num_stations_to_test -= 1
    
    s = random.randrange(values["start"], values["end"])
    s_offset = s-values["start"] # may get snode in the "middle" of the original file
    
    # create random timesteps to check 
    for rndt in range(num_timesteps_to_test):
      y = random.randrange(2004, 2017)
      m = random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      t_offset = t - int(date2num(datetime(y,m,1,0,0), units=u, calendar=c)) + tSTART # 0-744, because cyclical every month folder
      
      # For the .mat files above, the time index is specified in the name of the key-string, so getting the data is straightforward.
      # But to get the timestep (t) in the .spc file we need to find the t_offset, which is relative to the start time of the .spc file. 
      
      dateprint = num2date(t+tSTART, units=u, calendar=c)
      
      sfilepath = data_location+"/"+str(y)+"/"+f"{m:02d}"+"/results/"+station_name+".spc"
      try: spcdata = swan.loadSpc(sfilepath, monthOnly=m)["spectra"]
      except: print(f"couldnt read {sfilepath}"); continue
      
      local_snodes = spcdata[t_offset, s_offset, :, :] # time, nodes, frequency, direction
      rmote_snodes = swan["spc", "spectra", s, t][0][0] # otherwise we get [[[ data ]]]
      
      try:
        np.testing.assert_array_equal(local_snodes, rmote_snodes)
        print(f"snode {s} (offset={s_offset}) - {station_name} at {dateprint}. local file shape={spcdata.shape}  t={t} (offset={t_offset})  OK")
      except AssertionError as ae:
        print(f"snode {s} (offset={s_offset}) - {station_name} at {dateprint}. local file shape={spcdata.shape}  t={t} (offset={t_offset}) does NOT match bucket data")
      
  
if __name__ == "__main__":
  test_uvicswan_mat()
  test_uvicswan_spc()
