from datetime import datetime
from netCDF4 import num2date, date2num
from netcdfswan import NetCDFSWAN
import random
import numpy as np
from scipy.io import loadmat
import json

# localOnly should be set to False
# this will need the location of the swan data on the server 

data_location = ??? # wherever your original data is (relative to test_netcdfswan.py or use absolute path)
  
# constants
u, c = "hours since 1970-01-01 00:00:00.0", "gregorian"
tSTART = int(date2num(datetime(2004,1,1,0,0), units=u, calendar=c))  
num_stations_to_test = 13
num_timesteps_to_test = 1
  
# create input json object
inputFile = "./BCSWANv5.json" # set localOnly to False
with open(inputFile, "r") as f:
  inputJson = json.load(f)
inputJson["showProgress"] = False
  
def test_uvicswan_mat():
  global num_timesteps_to_test
  
  # load from s3 bucket
  swan = NetCDFSWAN(inputJson)
  
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
    
    # create random timesteps to check
    for i in range(num_timesteps_to_test):
      y = 2004 #random.randrange(2004, 2017)
      m = 1 #random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      dateprint = num2date(t+tSTART, units=u, calendar=c)
      
      mfilepath = data_location + str(y)+"/"+f"{m:02d}"+"/results/"+mVAR+".mat"
      try: matdata = loadmat(mfilepath)
      except: print(f"couldnt read {mfilepath}"); continue
      
      key = mvar+"_"+str(y)+f"{m:02d}"+f"{d:02d}"+"_"+f"{h:02d}"+"0000"
      local_nodes = matdata[key][0]
      rmote_nodes = swan["s", var, t][0]
      
      #print(rmote_nodes)
      
      np.testing.assert_array_equal(local_nodes, rmote_nodes)
      print(f"{key} {dateprint} OK")


def test_uvicswan_spc():
  global num_stations_to_test, num_timesteps_to_test
  
  # load from s3 bucket
  swan = NetCDFSWAN(inputJson)
  stn_slices = json.loads(swan.info()["metadata"]["station_slices"]) # need from aggregate of meta data
  
  #print(stn_slices)
  
  # check random snodes instead of stations
  for stn_key, snode_indices in stn_slices.items():
    if num_stations_to_test <= 0: break
    num_stations_to_test -= 1

    filename, snStart, snEnd = stn_key+".spc", snode_indices[0], snode_indices[1]
    s = random.randrange(snStart, snEnd)
    #s = 0
    
    # create random timesteps to check 
    for rndt in range(num_timesteps_to_test):
      y = 2004 #random.randrange(2004, 2017)
      m = 1 #random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      t_offset = t - int(date2num(datetime(y,m,1,0,0), units=u, calendar=c)) + tSTART # 0-744, because cyclical every month folder
      # For the .mat files above, the time index is specified in the name of the key-string, 
      # so its straightforward to get the data for it.
      # But to get the timestep (t) in the .spc file we need to find the t_offset, 
      # which is relative to the start time of the .spc file. 
      
      dateprint = num2date(t+tSTART, units=u, calendar=c)
      sfilepath = data_location + str(y)+"/"+f"{m:02d}"+"/results/"+filename
      
      spcdata = swan.loadSpc(sfilepath)
      spc = spcdata["spectra"]
      
      print(f"spc.shape={spc.shape} t={t}  t_offset={t_offset}  snode index={s}")
      
      local_snodes = spc[t_offset, s, :, :] # time, nodes, frequency, direction
      rmote_snodes = swan["spc", "spectra", s, t][0][0] # otherwise we get [[[ data ]]]
      #rmote_snodes = swan.query({"variable":"spectra","time":t,"snode":snStart}) # doesn't work
      
      #print(f"rmote_snodes.shape={rmote_snodes.shape}   local_snodes.shape={local_snodes.shape}")
      
      try:
        np.testing.assert_array_equal(local_snodes, rmote_snodes)
        print(f"snode {s} ({stn_key}) {dateprint}  OK")
      except AssertionError as ae:
        print(f"snode {s} ({stn_key}) {dateprint} does NOT match bucket data:\n{ae}")
  
if __name__ == "__main__":
  #test_uvicswan_mat()
  test_uvicswan_spc()
