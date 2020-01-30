#!/usr/bin/env python3
import os
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from createnodes import NodeMap

# note: pytest will not show progress bar if used in pycharm!
TQDM = True
try:
  from tqdm import tqdm
except ModuleNotFoundError:
  TQDM = False
  
TEST_FILE_PATH = "nc_examples/test.nc"
# could store in a different file?
TEMPORAL_VARIABLES = {
  "u": {
    "units": "m/s", 
    "standard name": "u velocity", 
    "long name": "u current velocity"
  },
  "v": {
    "units": "m/s", 
    "standard name": "v velocity", 
    "long name": "v current velocity"
  },
  "d": {
    "units": "m", 
    "standard name": "depth", 
    "long name": "water depth, in m (not the bottom level!)"
  },
  "h": {
    "units": "m", 
    "standard name": "?Free surface elevation", 
    "long name": "?Free surface elevation, m"
  },
  "u10": {
    "units": "m/s", 
    "standard name": "u Wind velocity", 
    "long name": "u Wind velocity, m/s"
  },
  "v10": {
    "units": "m/s", 
    "standard name": "v Wind velocity", 
    "long name": "v Wind velocity, m/s"
  },
  "hs": {
    "units": "m", 
    "standard name": "significant wave height", 
    "long name": "significant wave height (in s)"
  },
  "tps": {
    "units": "s", 
    "standard name": "peak period", 
    "long name": "\'smoothed\' peak period (in s)"
  },
  "tmm10": {
    "units": "s", 
    "standard name": "mean absolute wave period", 
    "long name": "mean absolute wave period (in s)"
  },
  "tm01": {
    "units": "s", 
    "standard name": "mean absolute wave period", 
    "long name": "mean absolute wave period (in s)"
  },
  "tm02": {
    "units": "s", 
    "standard name": "mean absolute zero-crossing period", 
    "long name": "mean absolute zero-crossing period (in s)"
  },
  "pdir": {
    "units": "degrees", 
    "standard name": "peak wave direction", 
    "long name": "peak wave direction in degrees"
  },
  "dir": {
    "units": "", 
    "standard name": "mean wave direction", 
    "long name": "mean wave direction (Cartesian or Nautical convention)"
  },
  "dspr": {
    "units": "degrees", 
    "standard name": "directional wave spread", 
    "long name": "directional spreading of the waves (in degrees)"
  },
  "qp": {
    "units": "", 
    "standard name": 
    "peakedness of wave spectrum", "long name": "peakedness of the wave spectrum (dimensionless)"
  },
  "transp": {
    "units": "m3/s", 
    "standard name": "transport of energy", 
    "long name": "transport of energy (in W/m or m3/s)"
  },
  "pthsign": {
    "units": "m", 
    "standard name": "partition of significant wave height", 
    "long name": "user requests partition of the significant wave height (in m)"
  },
  "ptrtp": {
    "units": "s", 
    "standard name": "partition of relative peak period", 
    "long name": "user requests partition of the relative peak period (in s)"
  },
  "ptdspr": {
    "units": "degrees", 
    "standard name": "partition of spread direction", 
    "long name": "user requests partition of the directional spreading (in degrees)"
  },
  "ptwlen": {
    "units": "m", 
    "standard name": "partition of average wave length", 
    "long name": "user requests partition of the average wave length (in m)"
  },
  "ptdir": {
    "units": "degrees", 
    "standard name": "partition of peak wave direction", 
    "long name": "user requests partition of the peak wave direction in degrees"
  },
  "ptsteep": {
    "units": "", 
    "standard name": "partition of wave steepness", 
    "long name": "user requests partition of the wave steepness (dimensionless)"
  }
}


def createNetCDF(filepath, node_map):
  with Dataset(filepath, "w") as src_file:
    
    src_file.title = "File description"
    src_file.institution = "Specifies where the original data was produced"
    src_file.source = "The method of production of the original data"
    src_file.history = "Provides an audit trail for modifications to the original data"
    src_file.references = "Published or web-based references that describe the data or methods used to produce it"
    src_file.comment = "Miscellaneous information about the data or methods used to produce it"
  
    # Dimensions
    ntime = src_file.createDimension("ntime", None)
    nnode = src_file.createDimension("nnode", node_map.num_nodes)
    npe = src_file.createDimension("npe", 3)
    nelem = src_file.createDimension("nelem", len(node_map.elements))
    
    #nnode_c = src_file.createDimension("nnode_c", 2)
    #nspectra = src_file.createDimension("nspectra", 50)
    #nfreq = src_file.createDimension("nfreq", 33)
    #ndir = src_file.createDimension("ndir", 36)
  
    # ------ Static variables
    elem = src_file.createVariable("elem","i4",("nelem","npe",))
    elem.units = ""
    elem.standard_name = "element"
    elem.long_name = "element"
  
    lon = src_file.createVariable("lon","f8", ("nnode",))
    lon.units = "degrees_east"
    lon.standard_name = "longitude"
    lon.long_name = "longitude"
  
    lat = src_file.createVariable("lat", "f8", ("nnode",))
    lat.units = "degrees_north"
    lat.standard_name = "latitude"
    lat.long_name = "latitude"
    
    
    # for sorting latitudes and longitudes. The assumption is that it will make
    # searching for nodes in a specified region quicker. (still TODO)
    lonsort = src_file.createVariable("lonsort","i4", ("nnode",))
    lonsort.units = "node number"
    lonsort.standard_name = "node number"
    lonsort.long_name = "node number of longitude"
  
    latsort = src_file.createVariable("latsort", "i4", ("nnode",))
    latsort.units = "node number"
    latsort.standard_name = "node number"
    latsort.long_name = "node number of latitude"
    
    
    time = src_file.createVariable("time", "f8", ("ntime",))
    time.units = "hours since 0001-01-01 00:00:00.0" 
    time.calendar = "gregorian"
    time.standard_name = "Datetime"
    time.long_name = "Datetime"
  
    # 'depth' of sea floor doesn't change. node_map.bathymetry
    b = src_file.createVariable("b", "f8", ("nnode",))
    b.units = "m"
    b.standard_name = "Bathymetry"
    b.long_name = "Bathymetry, m (CGVD28)"

    # ------ Temporal variables (have similar structure)
    for k,v in TEMPORAL_VARIABLES.items():
      var = src_file.createVariable(k, "f8", ("ntime", "nnode",))
      var.units = v['units']
      var.standard_name = v['standard name']
      var.long_name = v['long name']
  
  
    """
    # spectra component
    specnodeindex = src_file.createVariable("specnodeindex", "i4", ("nspectra",))
    # TODO
    spectra = src_file.createVariable("spectra", "f4", ("ntime","nspectra","nfreq","ndir",))
    spectra.units = "TODO"
    spectra.standard_name = "TODO"
    spectra.long_name = "Spectra"
    """
  

def writeNetCDF(filepath, node_map):
  with Dataset(filepath, "r+") as src_file:
    # Dimensions
    
    ntime = len(src_file.dimensions['ntime'])
    nnode = len(src_file.dimensions['nnode'])
    npe = len(src_file.dimensions['npe'])
    nelem = len(src_file.dimensions['nelem'])
    #nspectra = len(src_file.dimensions["nspectra"])
    #nfreq = len(src_file.dimensions["nfreq"])
    #ndir = len(src_file.dimensions["ndir"])

    elem = src_file.variables['elem']
    lon = src_file.variables['lon']
    lat = src_file.variables['lat']
    lonsort = src_file.variables['lonsort']
    latsort = src_file.variables['latsort']
    time = src_file.variables['time']
    
    
    b = src_file.variables['b']
    u = src_file.variables['u'] # current vel
    v = src_file.variables['v'] # current vel
    d = src_file.variables['d']
    #h = src_file.variables['h'] #?
    u10 = src_file.variables['u10'] # wind vel
    v10 = src_file.variables['v10'] # wind vel
    hs = src_file.variables['hs']
    tps = src_file.variables['tps']
    tmm10 = src_file.variables['tmm10']
    tm01 = src_file.variables['tm01']
    tm02 = src_file.variables['tm02']
    pdir = src_file.variables['pdir']
    dir_ = src_file.variables['dir']
    dspr = src_file.variables['dspr']
    qp = src_file.variables['qp']
    transp = src_file.variables['transp'] #may have to split x/y
    pthsign = src_file.variables['pthsign']
    ptrtp = src_file.variables['ptrtp']
    ptwlen = src_file.variables['ptwlen']
    ptdir = src_file.variables['ptdir']
    ptdspr = src_file.variables['ptdspr']
    ptsteep = src_file.variables['ptsteep']
    
    #specnodeindex = src_file.variables['specnodeindex']
    #spectra = src_file.variables['spectra']

    # ------ Static
    # Should only load static data once.
    # Temporal data (below) added one at a time otherwise.
    # could be a whole separate netCDF file in the end. need to save space
    node_no = 0
    
    elem[:] = np.array([n for n in node_map.elements])
    #print(f"CHECK elements of node {node_no}:", elem[node_no])
    
    lon[:] = np.array([lo for lo in node_map.lon])
    #print(f"CHECK longitude of node {node_no}:", lon[node_no])
    
    lat[:] = np.array([la for la in node_map.lat])
    #print(f"CHECK latitude of node {node_no}:", node_map.lat[node_no])
    
    b[:] = np.array([bath for bath in node_map.bathymetry])
    #print(f"CHECK bathymetry of node {node_no}:", node_map.bathymetry[node_no])
    
    # testing
    lonsort[:] = np.array([n for n in node_map.lon_sort])
    #print(f"CHECK node of westernmost longitude:", lonsort[0], "at", lon[lonsort[0]])
    #print(f"CHECK node of easternmost longitude:", lonsort[-1], "at", lon[lonsort[-1]])
    
    latsort[:] = np.array([n for n in node_map.lat_sort])
    #print(f"CHECK node of southernmost latitude:", latsort[0], "at", lat[latsort[0]])
    #print(f"CHECK node of northernmost latitude:", latsort[-1], "at", lat[latsort[-1]])
    
    
    #specnodeindex[:] = np.arange(0, nspectra, dtype=np.int32)
    
    # ------ Temporal (extra dimension!)
    t = 0
    
    # extract date elements with node_map.timesteps
    # change t[9:11] (hours) if necessary 
    # do we want to enumerate with i,
    # or use dict keys instead? e.g. '20040101_020000': 2004-1-1 2:0:0
    for i, ts in enumerate(node_map.timesteps):
      date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
      time[i] = date2num(date, units=time.units, calendar=time.calendar)
    print(f"CHECK time step for 1 Jan 2004:", time[0])
    
    #print("node ordering:", node_map.node_order)
    
    # store time coord as dictionary? 
    # faster access instead of a loop or list
    # e.g. '20040101_020000': -56.5
    # date2num becomes 17557968.0 for 1 Jan 2004
    # comes from  node_map.timesteps
    
    #u[t] = np.arange(0, nnode, dtype=np.float64) * 0.1
    #v[t] = np.arange(0, nnode, dtype=np.float64) * 0.1
    #d[t] = np.arange(0, nnode, dtype=np.float64) * 0.1
    #h[t] = np.arange(0, nnode, dtype=np.float64) * 0.1
    
    if TQDM:
      print("\n------  HS   ------")
      for tt in tqdm(range(len(time))):
        hs[tt, :] = np.array([n for n in node_map.matfiles['HS'][node_map.timesteps[tt]]])
      """
      print("\n------ WINDX ------")
      for tt in tqdm(range(len(time))):
        u10[tt, :] = np.array([n for n in node_map.matfiles['WINDX'][node_map.timesteps[tt]]])
      print("\n------ WINDY ------")
      for tt in tqdm(range(len(time))):
        v10[tt, :] = np.array([n for n in node_map.matfiles['WINDY'][node_map.timesteps[tt]]])
      print("\n------  TPS  ------")
      for tt in tqdm(range(len(time))):
        tps[tt, :] = np.array([n for n in node_map.matfiles['TPS'][node_map.timesteps[tt]]])
      print("\n------  QP   ------")
      for tt in tqdm(range(len(time))):
        qp[tt, :] = np.array([n for n in node_map.matfiles['QP'][node_map.timesteps[tt]]])
      print("\n------ DSPR  ------")
      for tt in tqdm(range(len(time))):
        dspr[tt, :] = np.array([n for n in node_map.matfiles['DSPR'][node_map.timesteps[tt]]])
      """
    else:
      for tt in range(len(time)):
        u10[tt, :] = np.array([n for n in node_map.matfiles['WINDX'][node_map.timesteps[tt]]])
    
      for tt in range(len(time)):
        v10[tt, :] = np.array([n for n in node_map.matfiles['WINDY'][node_map.timesteps[tt]]]) 
      
    print()
    #print("CHECK wind vel x of all nodes at time 0:", u10[0, :])
    #print("CHECK wind vel y of all nodes at time 0:", v10[0, :])
    # 24 * 31 + 1 = 745 timesteps
    


def printNetCDF(filepath):
  with Dataset(filepath, "r") as src_file:
    ncdump(src_file)

def ncdump(nc_fid, verb=True):
  '''
  ncdump outputs dimensions, variables and their attribute information.
  The information is similar to that of NCAR's ncdump utility.
  ncdump requires a valid instance of Dataset.

  Parameters
  ----------
  nc_fid : netCDF4.Dataset
      A netCDF4 dateset object
  verb : Boolean
      whether or not nc_attrs, nc_dims, and nc_vars are printed

  Returns
  -------
  nc_attrs : list
      A Python list of the NetCDF file global attributes
  nc_dims : list
      A Python list of the NetCDF file dimensions
  nc_vars : list
      A Python list of the NetCDF file variables
  '''
  
  def print_ncattr(key):
    """
    Prints the NetCDF file attributes for a given key

    Parameters
    ----------
    key : unicode
        a valid netCDF4.Dataset.variables key
    """
    try:
      print("type: " + repr(nc_fid.variables[key].dtype))
      for ncattr in nc_fid.variables[key].ncattrs():
        print('%s:'.format(ncattr)+repr(nc_fid.variables[key].getncattr(ncattr)))
    except KeyError:
      print("WARNING: %s does not contain variable attributes" % key)
  
  # NetCDF global attributes
  nc_attrs = nc_fid.ncattrs()
  if verb:
    print("NetCDF Global Attributes:")
    for nc_attr in nc_attrs:
      print('%s: ' % nc_attr + repr(nc_fid.getncattr(nc_attr)))
  nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
  
  # Dimension shape information.
  if verb:
    print("\nNetCDF dimension information:")
    for dim in nc_dims:
      print("Name:" + dim)
      print("size:" + str(len(nc_fid.dimensions[dim])))
      print_ncattr(dim)
  
  # Variable information.
  nc_vars = [var for var in nc_fid.variables]  # list of nc variables
  if verb:
    print("\nNetCDF variable information:")
    for var in nc_vars:
      if var not in nc_dims:
        print("\nName:" + var)
        # print("dimensions:" + str(nc_fid.variables[var].dimensions))
        print(nc_fid.variables[var].dimensions)
        # 745 timesteps x 177945 nodes = 132569025
        print(nc_fid.variables[var].size) 
        
        print_ncattr(var)
  return nc_attrs, nc_dims, nc_vars


    
  


if __name__ == "__main__":
  createNetCDF(TEST_FILE_PATH)
  writeNetCDF(TEST_FILE_PATH)
  printNetCDF(TEST_FILE_PATH)

