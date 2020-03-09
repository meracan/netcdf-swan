#!/usr/bin/env python3
"""
    netcdfswan.py
        Creates and stores partitioned netcdf (.nc) files from SWAN data into a local cache and
    uploads them to an Amazon S3 bucket, using the s3-netcdf (tentative, will change) package.

    The directory path for the SWAN data should be written in 'input_master.json' ahead of time so the script knows
    where to look when it runs. Other path names are also stored in this json file. Ideally, the script should be placed
    beside the data folder. netcdfswan.py assumes a folder structure sorted temporally,
    with 12 month folders per year folder:

    `-- SWAN_DATA
        |
        |-- 2005
        |-- 2006
        |-- 2007
        :   |-- 01
            |-- 02
            |-- 03
            :   `-- results
                    |-- HS.mat
                    |-- WIND.mat
                    |-- line_n.spc
                    |-- QP.mat
                    :

                    :
                    |-- hotspots.spc
                    `-- TPS.mat
            :
            |-- 11
            `-- 12
        :
        |-- 2016
        |-- 2017
        |
        `-- Mesh
            |-- .ele
            |-- .bot
            `-- .node

    The Mesh folder holds latitudes, longitudes, and bathymetry of all the nodes, and other triangle mesh information.
    The Matlab (.mat) files in the 'results' folder of each month holds variable information for all of the nodes,
    as well as spectra information.

    Procedure:
    swan data -> node map -> .nc + .nca files (NetCDF2D: write) -> upload (NetCDF2D: S3 client)

    The Node Map class acts as a sort of hub to hold all of the static and meta data as nc files are being created,
    and provides some useful methods for plotting, etc. When the Node Map is initialized, it first loads the static data
    contained in the Mesh folder.

    Uploading:
    1. Creates a node map that first searches for the Mesh folder in the path 'data_folder' and
        stores its .ele, .bot, and .node data.

    2. loads .mat files in the 'results' folder for each year and month

    3. uploads data to cache and cloud using s3netcdf package

    Basic usage:
    create node map first, then upload the .mat files:

    nm = NodeMap()
    nm.upload_mats()


    ...

    TODO:
        input_master: change localOnly to False for actual run
        upload(): give stats, how much data was transferred/uploaded, which data/dates, etc.
            'verbose' option
        download(): need? based on specificity, including variable:
            $ netcdf-swan.py download <year> <month> <variable>
            useful for visual manipulation, triangle mesh grids, .nc creation, etc.
        option to convert to non-partitioned .nc file



"""
import json
import os, sys
import numpy as np
import re
import pprint as pp
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
from s3netcdf import NetCDF2D

TQDM = True  # <--- change
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    TQDM = False

# regex for possible matlab file names.
re_mat = re.compile(r"(WIND|HS|TPS|TMM10|TM01|TM02|PDIR|DIR|DSPR|QP|TRANSP)[.]mat")
re_spc = re.compile(r"()[.]spc")
re_mat_test = re.compile(r"(WIND_6hr|HS_6hr|HS_)[.]mat")  # separate test from actual?

re_ele = re.compile(r".+[.]ele")
re_bot = re.compile(r".+[.]bot")
re_node = re.compile(r".+[.]node")
re_x = re.compile(r"^.*_x_")
re_y = re.compile(r"^.*_y_")
re_upload = re.compile(r"^-?[uU](pload)?$")
re_download = re.compile(r"^-?[dD](ownload)?$")
re_static = re.compile(r"^-?[sS](tatic)?$")
re_mesh = re.compile(r"^-?[mM](esh)?$") #
re_edit = re.compile(r"^-?[eE](dit)?$")  # might not use, but we'll see
re_year = re.compile(r"^20[0-9][0-9]$")
re_month = re.compile(r"^(0?[1-9]|1[0-2])$")
# re_day = re.compile(r"^(0?[1-9]|[12][1-9]|3[01])$")
# re_hour = re.compile(r"^(0?[0-9]|1[0-9]|2[0-3])(00)?$")
# re_lat = re.compile(r"^[-]?lat$")
# re_lon = re.compile(r"^[-]?long?$")


class NodeMap:

    def __init__(self, reloading=False):
        """ 1. reads and stores .ele, .bot, and .node data in the Mesh folder
            2. creates timesteps based on the earliest and latest possible times in the dataset.
            3. creates the nca template using the mesh and time information
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


    def load_mat(self, filepath, mname):
        """
            Loads one mat file into self.matfile1 and self.matfile2 if there are both X and Y coordinates
            (e.g. WIND.mat), otherwise just stores into self.matfile1.

            The header, version and globals are popped for now. (might need later?)

            (self.timesteps used to be replaced each time)

            'np.squeeze()' is needed to get rid of nested arrays like 'array([[-22.5, -18.0, -36.7, ...]])'

            name stripping is to isolate the date:
                key <- 'abcde_20040101_050000'
                key[-15:] -> '20040101_050000'

            the start date is updated with the first time step of the month at the end
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

    # ===================== below not needed anymore?
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


    # ===========================


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
        - Pseudocode -

        determine starting month folder
        for each year:
            for each month:
                for each .mat (and .spc) file:
                    load .mat file into the Node Map (e.g. HS.mat)
                    use s3-netcdf package to partition the data and store in the cache
                    set reloading to false to prevent mesh data from loading again (only need to do once)
                write latest timestep to 'startdate.txt'

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




    def create_nca_input(self):
        """
        Creates an 'nca' structure for the master input template (json) using the swan mesh data.
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


    def upload(self):
        """
            A NetCDF2D object is created with the finished template and the contents from the node map
                are written to the NetCDF2D object one timestep at a time.
            if 'reloading' is true, that means we need to re-load all static data as well as the time steps.

            upload_to_cache then empties both self.matfile1 and self.matfile2 --> {} for the next upload
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




    def download(self):
        """
        TODO
        Create NetCDF2D object like the upload method, but use it to read from cache instead.
        Can then use tri_plot to visualize the data.

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
        changed load_mat to only get one matfile at a time, then uploading.
        Python only has a limited memory capacity, and since each mat file is just over half a gigabyte,
        the matfiles = {} dictionary will not be able to hold all of the information:

            0 - mesh data only needs to be read once, when calling NodeMap() for the first time.

            upload_mats():
              for year in years:
                for month in months:
                  1 - load ONE mat file in the NodeMap (e.g. HS.mat)
                  3 - upload mat file to the cloud
                  4 - delete mat file from nodemap
                  5 - start from step 1 with NEW mat file (e.g. WIND.mat)

        After a mat file is loaded into the node map, meta data is put into an input template for the
        master nc file (.nca). The .nca file template is then updated with the timesteps that
        were extracted from the .mat data.


        should have parameter to indicate where to continue from if upload is interrupted
        In this case the folders are assumed to be extracted from in chronological order.

        might make more sense if folder path doesn't need to be specified....

        if one of the variables fails to upload, whole month should be re-uploaded again (not too inconvenient)
        basically, can start from any year or month of a year and upload from those onward


        if upload does get interrupted, should be able to detect where it left off automatically
            based on timestep information in the nca file
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

        $ python3 netcdfswan.py upload mesh  # deprecated


        # putting in dates implies the timesteps need to be loaded again
        $ python3 netcdfswan.py upload 2004 01
        $ python3 netcdfswan.py upload 2004
        $ python3 netcdfswan.py upload

        $ python3 netcdfswan.py download 2007 04 hs # dont worry about until after

    """
    # get rid of parentheses and script name
    args = args[0]
    args.pop(0)

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
    #print(args_)
    main(args_)