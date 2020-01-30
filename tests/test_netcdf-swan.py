#!/usr/bin/env python3
import sys, os
import re
from createnodes import NodeMap
from createnetcdf import createNetCDF, writeNetCDF, printNetCDF

TEST_FILE_PATH = "../output/test.nc"
MESH_FOLDER = "../data/Mesh"
XYZ = "../data/XYZ.mat"
TEMP_2004_01_RESULTS = "../data/2004/01/results"
TEMP_2004_01_RESULTS_WIND = "../data/2004/01/results/WIND.mat"

def test_main():
  """ swan data -> node map -> netcdf -> s3 object -> deploy
      may shortcut some steps above (swan data -> netcdf)

      1. Searches for XYZ and mesh files in the path directory specified by the user (no arguments or '.' 
        means search current directory) and creates a temporary node map.
      Arguments after that are date and time coordinates, to search a specific month of a particular year 
        for the 'results' folder, e.g.   $ python3 netcdf-swan.py . 2004 01
        will find the 'results' data folder in January of 2004, in the current directory.
      If no date-time is specified, should have a loop to get ALL results folders. (TODO)
      
      2. Reads all .mat files in the results folder and stores them as node data.
      3. Creates and writes a netcdf file (.nc) with the node data
  
      4. read node data from nc file (for testing)
      5. plot ?
      
  
  """
  xyz, mesh = XYZ, MESH_FOLDER
  results = TEMP_2004_01_RESULTS # temporary
  
  args = sys.argv
  if len(args) > 1:
    filepath = sys.argv[1]
    if os.path.isdir(filepath):
      for root, dirs, files in os.walk(filepath):
        if "XYZ.mat" in files:
          xyz = os.path.join(root, "XYZ.mat")
        if "Mesh" in dirs:
          mesh = os.path.join(root, "Mesh")
      if len(args) > 3:
        results = os.path.join(args[1], args[2], args[3], "results")
        #print("results path:", results)
  
  print(f"\n*** Loading node data from \'{results}\' ***")
  nm = NodeMap()
  nm.load_mesh(mesh, xyz)
  nm.load_mat(results)
  nm.sort_coords()
  
  
  
  print(f"\n*** Creating NetCDF file as \'{TEST_FILE_PATH}\' ***")
  createNetCDF(TEST_FILE_PATH, nm)
  writeNetCDF(TEST_FILE_PATH, nm)
  
  
  print(f"\n*** (testing) Reading created NetCDF file ***")
  nm_from_nc = NodeMap()
  nm_from_nc.load_nc(TEST_FILE_PATH)
  
  box = nm_from_nc.get_node_area_delta(-123.3656, 48.4284, 1.0)
  box2 = nm_from_nc.get_node_area([-124.3656, 47.4284, -123.3656, 48.4284])

  print(box[:10], "...")
  print("number of nodes around victoria (1 degree box):", len(box2))

  
  # pyplot image
  # nm.tri_plot()
  # nm.update_plot()
  
  print("\n*** Finished ***")
  
  


if __name__ == "__main__":
  test_main()
  