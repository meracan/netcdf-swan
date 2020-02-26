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

    Uploading:
    1. Creates a node map that first searches for the Mesh folder in the path 'data_folder' and
        stores its .ele, .bot, and .node data.

    2. loads mat files in the 'results' folder for each given year and month if specified.

    3. uploads data to cache and cloud using s3netcdf package

    Downloading:
    ???

"""
import json
import os, sys
import numpy as np
import re
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
from s3netcdf import NetCDF2D
# from createnetcdf import create_nca_input

TQDM = True  # <--- change
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
re_upload = re.compile(r"-?[uU](pload)?")
re_download = re.compile(r"-?[dD](ownload)?")
re_edit = re.compile(r"-?[eE](dit)?")  # might not use, but we'll see
re_year = re.compile(r"^20[0-9][0-9]$")
re_month = re.compile(r"^(0?[1-9]|1[0-2])$")
# re_day = re.compile(r"^(0?[1-9]|[12][1-9]|3[01])$")
# re_hour = re.compile(r"^(0?[0-9]|1[0-9]|2[0-3])(00)?$")
# re_lat = re.compile(r"^[-]?lat$")
# re_lon = re.compile(r"^[-]?long?$")


class NodeMap:

    def __init__(self):
        """
            only plot a subset of all the node data, which should be loaded first (change later)
            dynamic plot, visual rep. of the nodes
            needs latitude and longitude (bounding box)

            When initialized, the NodeMap automatically searches for the Mesh folder in the
                path 'data_folder' and stores its .ele, .bot, and .node data.
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

        with open("../filepath_names.json") as fpnames:
            names = json.load(fpnames)
            self.mesh_folder = names['mesh_folder']
            self.data_folder = names['data_folder']  # ./BCSWANv5, but uses test_variable_names anyway
        with open("../test_variable_names.json") as vnames: # !!! remove 'test' when finished
            self.temp_var_names = json.load(vnames)
        with open("../read_master.json") as rt:
            self.read_template = json.load(rt)
        with open("../input_master.json") as im:
            self.master_input = json.load(im)

        # mesh data and nca input are loaded and created automatically
        self.load_mesh(self.mesh_folder) # will put num_nodes and num_elements
        self.create_nca_input()


    def load_mesh(self, meshfolder):
        """ loads all mesh data and stores into class object
            assumes mesh folder exists
            (should need to do only once--)
            * UPDATE
            should do automatically when NodeMap() is called
        """
        mesh_node = ""
        mesh_ele = ""
        mesh_bot = ""
        for name in os.listdir(meshfolder): # necessary?
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
            otherwise the filename should look like '[...]/results/WIND.mat' or something (not working yet)

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
        print(filename)
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


        # update timesteps in input template
        self.master_input["nca"]["dimensions"]["ntime"] = len(self.timesteps)


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
    # must be a clean, empty node map to use the below methods (not used anymore?)


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
            for k, v in self.temp_var_names.items():
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



    def load_mats(self, year="", month=""):
        """
        After mat files are loaded into the node map, meta data is put into an input template for the
        master nc file (.nca). The .nca file template is then updated with the timesteps that
        were extracted from the .mat data.

        """
        print(f"*** loading mat data")
        data_folder = self.data_folder
        if year:
            if month:  # .mat files for one month
                results = os.path.join(data_folder, year, month, "results")
            else:  # .mat files for all months of one year
                results = os.path.join(data_folder, year)

            self.load_mat(results)
        else:  # .mat files for ALL months from ALL years
            for year in os.listdir(data_folder):
                if year != 'Mesh' and not year.startswith('.') and os.path.isdir(year):  # avoid '.DS_Store' file
                    for month in os.listdir(data_folder + "/" + year):
                        results = os.path.join(data_folder, year, month, "results")
                        self.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

        #print("*** creating nca template and adding timesteps")
        #self.master_input["nca"]["dimensions"]["ntime"] = len(self.timesteps)


    def create_nca_input(self):  # to eliminate 'createnetcdf.py', probably only need the one .py script
        """
        Creates an 'nca' structure for the master input template (json) using swan mesh data.
        No variable data is loaded yet.

        Information stored in another .json file
        """

        master_input = self.master_input

        metadata = dict(
            title="File description",
            institution="Specifies where the original data was produced",
            source="The method of production of the original data",
            history="Provides an audit trail for modifications to the original data",
            references="Published or web-based references that describe the data or methods used to produce it",
            comment="Miscellaneous information about the data or methods used to produce it"
        )
        master_input["metadata"] = metadata

        # no timesteps (ntime) in the static nc file, will be updated as more temporal nc files are created
        dimensions = dict(
            npe=3,
            nnode=self.num_nodes,
            ntime=0,
            nelem=self.num_elements
        )
        master_input["nca"]["dimensions"] = dimensions

        variables = dict(
            lat=dict(type="float64", dimensions=["nnode"], units="degrees_north", standard_name="latitude",
                     long_name="latitude"),
            lon=dict(type="float64", dimensions=["nnode"], units="degrees_east", standard_name="longitude",
                     long_name="longitude"),
            elem=dict(type="int32", dimensions=["nelem", "npe"], units="", standard_name="elements",
                      long_name="elements"),
            time=dict(type="float64", dimensions=["ntime"], units="hours since 1970-01-01 00:00:00.0",
                      calendar="gregorian",
                      standard_name="Datetime", long_name="Datetime"),
            b=dict(type="float32", dimensions=["nnode"], units="m", standard_name="Bathymetry",
                   long_name="Bathymetry, m (CGVD28)")
        )
        master_input["nca"]["variables"] = variables

        # temporal data (check for master)
        variables2 = {}
        for k, v in self.temp_var_names.items():
            var = dict(type="float32", units=v["units"], standard_name=v["standard name"], long_name=v["long name"])
            variables2.update({k: var})

        groups = dict(
            s=dict(dimensions=["ntime", "nnode"], variables=variables2),
        )
        master_input["nca"]["groups"] = groups

        self.master_input = master_input
        #return master_input


    def upload(self):
        """
        A NetCDF2D object is created with the finished template and the contents from the node map
            are written to the NetCDF2D object one timestep at a time.
        """
        print("- Upload -")
        print("*** initializing NetCDF2D object")
        netcdf2d = NetCDF2D(self.master_input)

        print("*** loading nca data into NetCDF2D object")
        for kvar, val in self.temp_var_names.items(): #TEMPORAL_VARIABLES.items():
            MAT_val = val['matfile name']
            if MAT_val in self.matfiles:
                for kv, n in self.matfiles[MAT_val].items():
                    ts = kv[-15:]
                    date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                    i = self.timesteps.index(date)
                    try:
                        netcdf2d["s", kvar, i] = n
                    except ValueError:
                        raise

        print("*** NetCDF2d created and shipped")

        # nm.upload() # uploading to s3 using NodeMap() instead of the other script
        # will also use this for downloading as well, so of course need documentation to help whoever will be using it


    def download(self):
        """
        Create NetCDF2D object like the upload method, but use it to read from cache instead.
        Can then use tri_plot to visualize the data
        """
        print("- Download -")
        print("*** initializing NetCDF2D object")
        netcdf2d = NetCDF2D(self.read_template)


    def print_netcdf(self, filepath, verb=True):
        with Dataset(filepath, "r") as src_file:
            self.ncdump(src_file, verb)

    def ncdump(nc_fid, verb=True):
        """
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
        """

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
                    print('%s:'.format(ncattr) + repr(nc_fid.variables[key].getncattr(ncattr)))
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


def main(*args):
    """
        (User parameters)

        If a year and month are specified, will search the results folder just in that month:
            $ netcdf-swan.py upload <data folder path> <year> <month>

        Otherwise, will look in the results folder for each month of the year if only that year is specified:
            $ netcdf-swan.py upload <data folder path> <year>

        Otherwise if no path is specified, ALL results folders will be searched (default: no arguments or '.' means
            search current directory):
            $ netcdf-swan.py upload <data folder path>
                OR
            $ netcdf-swan.py upload

        e.g. To find the 'results' data folder in January of 2004 in the current directory:
            $ python3 netcdf-swan.py . 2004 01

        :param args:
            <upload/download> <folder path> <year> <month> ...
        :return:
    """

    args = args[0]  # get rid of parentheses

    i = 0
    uploading = re.match(re_upload, args[i]) is not None
    downloading = re.match(re_download, args[i]) is not None

    get_year, get_month, get_day, get_hour = "", "", "", ""
    i += 1
    if not re.match(re_year, args[i]) and not re.match(re_month, args[i]):  # if data folder path
        if os.path.isdir(args[i]):
            data_folder = args[i]
        i += 1
    if len(args) > i and re.match(re_year, args[i]):
        get_year = args[i]
        i += 1
    if len(args) > i and re.match(re_month, args[i]):
        get_month = args[i]
        i += 1

    if uploading:
        nm = NodeMap()
        nm.load_mats(get_year, get_month) # get_var?
        nm.upload()

    elif downloading:
        nm = NodeMap()
        nm.download()


    print("*** finished")


if __name__ == '__main__':
    args_ = sys.argv
    main(args_[1:])