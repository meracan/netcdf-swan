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
num_stations_to_test = 5
num_timesteps_to_test = 7

def test_uvicswan_mat():
  global num_timesteps_to_test
  
  # create input json object (needs nca object)
  inputFile = "./BCSWANv5.json"
  with open(inputFile, "r") as f:
    inputJson = json.load(f)
  inputJson["showProgress"] = False
  
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
      y = random.randrange(2004, 2017)
      m = random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      dateprint = num2date(t+tSTART, units=u, calendar=c)
      
      mfilepath = data_location + str(y)+"/"+f"{m:02d}"+"/results/"+mVAR+".mat"
      try: matdata = loadmat(mfilepath)
      except: print(f"couldnt read {mfilepath}"); continue
      
      key = mvar+"_"+str(y)+f"{m:02d}"+f"{d:02d}"+"_"+f"{h:02d}"+"0000"
      local_nodes = matdata[key][0]
      rmote_nodes = swan["s", var, t]
      
      np.testing.assert_array_equal(local_nodes, rmote_nodes)
      print(f"{key} {dateprint} OK")


def test_uvicswan_spc():
  global num_stations_to_test, num_timesteps_to_test
  
  # create input json object
  inputFile = "./BCSWANv5.json" # set localOnly to False
  with open(inputFile, "r") as f:
    inputJson = json.load(f)
  inputJson["showProgress"] = False

  # load from s3 bucket
  swan = NetCDFSWAN(inputJson)
  stns = json.loads(swan.info()["metadata"]["stations"])
  
  for s in stns.items():
    if num_stations_to_test <= 0: break
    num_stations_to_test -= 1

    filename, i = s[0]+".spc", s[1]

    # create random timesteps to check 
    for rndt in range(num_timesteps_to_test):
      y = random.randrange(2004, 2017)
      m = random.randrange(1, 13)
      d = random.randrange(1, 29)
      h = random.randrange(0, 24)
      t = int(date2num(datetime(y,m,d,h,0), units=u, calendar=c)) - tSTART
      t_offset = t - int(date2num(datetime(y,m,1,0,0), units=u, calendar=c)) + tSTART
      dateprint = num2date(t+tSTART, units=u, calendar=c)
      mfilepath = data_location + str(y)+"/"+f"{m:02d}"+"/results/"+filename
      
      spcdata = swan.loadSpc(mfilepath)
      spc = spcdata["spectra"]
      
      # ------- just use first node in station ---------
      data_local = spc[t_offset, 0, :, :]
      data_rmote = swan.query({"group":"spc","variable":"spectra","station":i,"time":t,"snode":0})
      try:
        np.testing.assert_array_equal(data_local, data_rmote)
        print(f"station {s[0]} {dateprint}  OK")
      except AssertionError as ae:
        print(f"station {s[0]} {dateprint} does NOT match bucket data:\n{ae}")
  
if __name__ == "__main__":
  test_uvicswan_mat()
  test_uvicswan_spc()
