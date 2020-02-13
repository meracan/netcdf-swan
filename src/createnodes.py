#!/usr/bin/env python3
"""
    used before creating and loading to netcdf file, as well as for reading from
    a netcdf file to be used for plotting, etc.


    - read from swan files into memory      <-- this step
    - store node info as netCDF file
    - turn into s3 bucket object and deploy

      .node file in the mesh will have lat and long information
    load Mesh files -- 
      ele for triangles
      node for lat and long (as well as type)
  
    start with small data for now
  
    import h5py to resolve mat files not loading properly (Error):
    https://stackoverflow.com/questions/17316880/reading-v-7-3-mat-file-in-python

    if path is not specified, assumes data is in a folder called 'data'

    1 netcdf for each timestep.


"""
import boto3
import json
import os, sys
import pandas as pd
import numpy as np
import re
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime

TQDM = True
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    TQDM = False

# regex for possible matlab file names.
re_mat = re.compile(r"(WIND|HS|TPS|TMM10|TM01|TM02|PDIR|DIR|DSPR|QP|TRANSP)[.]mat")
re_mat_test = re.compile(r"(WIND_6hr|HS_6hr)[.]mat")  # separate test from actual?

re_ele = re.compile(r".+[.]ele")
re_bot = re.compile(r".+[.]bot")
re_node = re.compile(r".+[.]node")
re_x = re.compile(r"^.*_x_")
re_y = re.compile(r"^.*_y_")
    
# temporary files/directories. will need to change.
# for testing this .py as executable
with open("../filepath_names.json") as fpnames:
    names = json.load(fpnames)
    SWAN_FOLDER = names['data_folder']
    MESH_FOLDER = names['mesh_folder']
    NC_FILE_MESH = names['output_nc_mesh']
    NC_FILE = names['output_nc']

with open("../variable_names.json") as vnames:
    TEMPORAL_VARIABLES = json.load(vnames)


class NodeMap:

    def __init__(self):
        """ only plot a subset of all the node data, which should be loaded first (change later)
            dynamic plot, visual rep. of the nodes
            needs latitude and longitude (bounding box)
        """
        self.num_nodes = 0
        self.num_elements = 0
        self.node_order = None
        self.lon, self.lat = None, None
        self.lon_sort, self.lat_sort = [], []
        self.boundary_markers = None
        self.matfiles = {}
        self.elements = None
        self.timesteps = []
        self.bathymetry = []


    def load_mesh(self, meshfolder):
        """ loads all mesh data and stores into class object
            assumes mesh folder exists
            (should need to do only once--)
        """
        mesh_node = ""
        mesh_ele = ""
        mesh_bot = ""
        for name in os.listdir(meshfolder):
          if re.match(re_node, name):
            mesh_node = name
          elif re.match(re_ele, name):
            mesh_ele = name
          elif re.match(re_bot, name):
            mesh_bot = name

        self.num_nodes, self.node_order, self.lon, self.lat, self.boundary_markers = self.swan_noderead(meshfolder+"/"+mesh_node)
        self.elements = self.swan_eleread(meshfolder+"/"+mesh_ele)
        self.num_elements = len(self.elements)
        self.bathymetry = self.swan_botread(meshfolder+"/"+mesh_bot)

    # For the swan methods below, "-1" turns MATLAB's 1 indexing into Python's 0 indexing
    def swan_noderead(self, nodefile):
        f = open(nodefile, 'r')
        V = np.loadtxt(f)
        f.close()
        num_nodes = int(V[0,0])
        node_order = V[1:, 0].astype(int) - 1
        lon, lat, =  V[1:,1], V[1:,2]
        boundary_marker = V[1:,3].astype(int)
        return num_nodes, node_order, lon, lat, boundary_marker

    def swan_eleread(self, elefile):
        f = open(elefile,'r')
        V = np.loadtxt(f,skiprows=1)
        f.close()
        elements = V[:,1:] - 1
        elements = elements.astype(int)
        return elements

    def swan_botread(self, botfile):
        f = open(botfile,'r')
        V = np.loadtxt(f)
        f.close()
        z = V
        return z


    def load_mat(self, filename):
        """ searches the entire folder if month folder is given,
            otherwise the filename should look like '[...]/results/WIND.mat' or something

            header, version and globals are popped for now, might need later?

            need  'np.squeeze()'  to get rid of nested arrays like 'array([[-22.5, -18.0, -36.7, ...]])'

            matfile has either one dimension or two (x and y coordinates) for each time step.

            name stripping:
                key == 'abcde_20040101_050000'
                key[-15:] == '20040101_050000'

            Hsig_20040101_050000 is one column, which should be turned into one .nc file
        """
        folder = ""
        mfiles = [filename]
        timestep_keys = []
        if os.path.isdir(filename):
            folder = filename+"/"
            mfiles = [name for name in os.listdir(filename)]  # e.g. ['HS.mat'] or ['HS.mat', 'WIND.mat', ... ]
        for mfilename in mfiles:
            if re.match(re_mat, mfilename) or re.match(re_mat_test, mfilename):  # test included in here for now
                print(f" reading {mfilename}...")
                matfile_dict = {}
                matfile_dictx = {}
                matfile_dicty = {}
                try:
                    matfile = loadmat(folder+mfilename)
                    matfile.pop('__header__')
                    matfile.pop('__version__')
                    matfile.pop('__globals__')

                    # k: 'Hsig_20040101_030000' or 'Windx_20040106_120000', etc
                    # v: 'array([-22.5, -18.0, -36.7, ...])' for all 177495 nodes
                    for k, v in matfile.items():
                        v_sq = np.squeeze(v)
                        ts = k[-15:]
                        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                        timestep_keys.append(date)

                        if re.match(re_x, k):
                            matfile_dictx.update({k: v_sq})
                        elif re.match(re_y, k):
                            matfile_dicty.update({k: v_sq})
                        else:
                            matfile_dict.update({k: v_sq})

                    if len(matfile_dict) == 0:
                        self.matfiles[mfilename[:-4]+"X"] = matfile_dictx
                        self.matfiles[mfilename[:-4]+"Y"] = matfile_dicty
                    else:
                        self.matfiles[mfilename[:-4]] = matfile_dict

                    ttt = self.timesteps
                    kkk = timestep_keys
                    ttt.extend(kkk)
                    ttt = sorted(list(set(ttt)))
                    self.timesteps = ttt

                except NotImplementedError as e:
                  print(f"{e}!")


    def sort_coords(self):
        """ should be after getting coordinates. check if tqdm is installed
        """
        lon_i_sort = np.argsort(self.lon[:])
        lat_i_sort = np.argsort(self.lat[:])

        if TQDM:
            print("\nsorting longitudes...")
            for i in tqdm(range(len(lon_i_sort))):
                self.lon_sort.append( lon_i_sort[i] )

            print("\nsorting latitudes...")
            for i in tqdm(range(len(lat_i_sort))):
                self.lat_sort.append( lat_i_sort[i] )
        else:
            print("\nsorting latitudes and longitudes...")
            for i in range(len(lon_i_sort)):
                self.lon_sort.append( lon_i_sort[i] )
            for i in range(len(lat_i_sort)):
                self.lat_sort.append( lat_i_sort[i] )



    # ====================== from .nc
    # must be a clean, empty node map to use the below methods


    def load_nc_static(self, filename):
        """
            reads static data variables from nc file
            ? self.boundary_markers = None
        """
        with Dataset(filename, "r", format="NETCDF4") as nc:
            time = nc.variables['time']
            lon, lat = nc.variables['lon'], nc.variables['lat']
            lonsort, latsort = nc.variables['lonsort'], nc.variables['latsort']

            self.num_nodes = len(lon)
            # turn into lists to speed up algorithm (doesnt run properly otherwise)
            self.lon, self.lat = list(lon[:]), list(lat[:])
            self.lon_sort, self.lat_sort = list(lonsort[:]), list(latsort[:])

            self.bathymetry = nc.variables['b'][:]
            self.timesteps = nc.variables['time'][:]
            self.ready_timesteps = True
            self.elements = nc.variables['elem'][:]

        #print(self.timesteps)

    def load_nc_temporal(self, filename):
        """
            for each temporal data variable in nc file, store into self.matfiles dictionary.
            could specify range of time steps? will need this functionality after grabbing from s3
        """
        with Dataset(filename, "r", format="NETCDF4") as nc:
            for k, v in TEMPORAL_VARIABLES.items():
                var = nc.variables[k]
                mname = v["matfile name"]
                try:
                    self.matfiles[mname] = nc.variables[var][:]
                except KeyError as e:
                    cause = e.args[0]
                    print(cause)

        #print(self.timesteps)



    def node_submap_area(self, min_lat, max_lat, min_long, max_long, start_time, end_time):
        """ should have read a node map already
            * * probably as lambda, after deploying to bucket
        """
        self.min_lat = min_lat
        self.max_lat = max_lat
        self.min_long = min_long
        self.max_long = max_long
        self.start_time = start_time
        self.end_time = end_time


    def tri_plot(self, x, y, z=None, node_area=None):
        """
            optional method, just to see
            should be able to change self.depth to other z attributes.
            if node area has been specified, use those coordinates instead
        """
        if z is None:
          z = self.bathymetry

        tg = Triangulation(x, y, self.elements) # will need to change elements if only a specified region
        plt.figure() # figsize=(10,10))
        plt.gca().set_aspect('equal')
        plt.tripcolor(tg, z, edgecolors='k', vmin=min(z), vmax=max(z))
        plt.colorbar()
        plt.title('Contour plot of user-specified triangulation')
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.show()

    def update_plot(self):
        """ dynamic, need some sort of signal or command to trigger this
        """
        print("<update not implemented yet!>")


    def get_node_area(self, box=[0,0,0,0]):
        """ Returns a list of nodes within a boxed region: [west, south, east, north]
        """
        lon_box = [l for i,l in enumerate(self.lon_sort) if (box[0] <= self.lon[i] and self.lon[i] <= box[2])]
        lat_box = [l for i,l in enumerate(self.lat_sort) if (box[1] <= self.lat[i] and self.lat[i] <= box[3])]
        box = list(set(lon_box) & set(lat_box))
        return box

    def get_node_area_delta(self, longitude, latitude, delta):
        """ Returns a list of nodes within a region. The point is at a longitude
            and latitude, with each side at +/- delta degrees from that point
        """
        lon_box = [l for i,l in enumerate(self.lon_sort) if (longitude-delta <= self.lon[i] and self.lon[i] <= longitude+delta)]
        lat_box = [l for i,l in enumerate(self.lat_sort) if (latitude-delta <= self.lat[i] and self.lat[i] <= latitude+delta)]
        box = list(set(lon_box) & set(lat_box))
        return box


def main():
  """ 
  """
  
  # initialize node map
  #nm = node_map()
  
  # mesh folder and XYZ.mat file
  #nm.load_mesh(MESH_FOLDER, XYZ)
  
  # results folder should have all the other .mat files
  #nm.load_matfiles(TEMP_FOLDER_2004_01_RESULTS)
  
  
  #nm.get_timesteps() #??
  
  # show what it looks like (for now)
  # nm.tri_plot()
  
  # just for consistency
  # nm.update_plot()
  
  
  #print("success")
  
  

if __name__ == '__main__':
  main()
  
  #mesh_folder = sys.argv[1]
  #mat_files = sys.argv[2]