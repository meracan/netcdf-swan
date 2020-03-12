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
    nm.upload_files()


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
import os, sys, time
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
re_spc = re.compile(r"(c_dixon|c_eliz|campbell|e_dell|hotspots|line_n|line_s|line_w|m_nomad|"
                    r"n_hecat|ne_isle|neah|p_renf|perouse|s_hecat|s_morsb|s_nomad|sombrio|"
                    r"sooke|tarbotn|tillamk|tofino|w_dixon|w_morsb|w_otter|w_washn)[.]spc")

re_spcdate = re.compile(r"^20[0-9]{6}[.][0-9]{6}$")
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

    def __init__(self):
        """ 1. reads and stores .ele, .bot, and .node data in the Mesh folder
            2. creates timesteps based on the earliest and latest possible times in the dataset.
            3. creates the nca template using the mesh and time information
        """
        self.num_nodes = 0
        self.num_elements = 0
        self.num_timesteps = 0
        self.node_order = None
        self.lon, self.lat = None, None
        self.lon_sort, self.lat_sort = [], []
        self.boundary_markers = None
        self.matfile1 = {}
        self.matfile2 = {}
        self.spcfile = {}
        self.afreqs = []
        self.ndir = []
        self.elements = None
        self.timesteps = []
        self.start_date = ""
        self.start_year = ""
        self.start_month = ""
        self.bathymetry = []

        with open("../jsons/variable_names.json") as vnames:
            self.temp_var_names = json.load(vnames)
        with open("../jsons/read_master.json") as rt:
            self.read_template = json.load(rt)
        with open("../jsons/input_master.json") as im:
            self.master_input = json.load(im)
            self.mesh_folder = self.master_input["file paths"]["mesh folder"]
            self.data_folder = self.master_input["file paths"]["data folder"]
            self.start_year = self.master_input["start year"]
            self.start_month = self.master_input["start month"]

        # below could be optional if this class is used more generally/universally
        self.load_mesh()
        self.load_timesteps()  # user specified (?) if interruption, will need to read from s3 first.
        self.create_nca_input()  # timesteps can go in here now that all of them are 'loaded'


    def load_mesh(self):
        # For the swan methods below, "-1" turns MATLAB's 1 indexing into Python's 0 indexing
        def swan_noderead(self, nodefile):
            with open(nodefile, 'r') as f:
                V = np.loadtxt(f)
            num_nodes = int(V[0, 0])
            node_order = V[1:, 0].astype(int) - 1
            lon, lat, = V[1:, 1], V[1:, 2]
            boundary_marker = V[1:, 3].astype(int)
            return num_nodes, node_order, lon, lat, boundary_marker

        def swan_eleread(self, elefile):
            with open(elefile, 'r') as f:
                V = np.loadtxt(f, skiprows=1)
            elements = V[:, 1:] - 1
            elements = elements.astype(int)
            return elements

        def swan_botread(self, botfile):
            with open(botfile, 'r') as f:
                V = np.loadtxt(f)
            z = V
            return z

        mesh_node = ""
        mesh_ele = ""
        mesh_bot = ""
        for name in os.listdir(self.mesh_folder):
          if re.match(re_node, name):
            mesh_node = name
          elif re.match(re_ele, name):
            mesh_ele = name
          elif re.match(re_bot, name):
            mesh_bot = name

        self.num_nodes, self.node_order, self.lon, self.lat, self.boundary_markers = swan_noderead(self.mesh_folder+"/"+mesh_node)
        self.elements = swan_eleread(self.mesh_folder+"/"+mesh_ele)
        self.num_elements = len(self.elements)
        self.bathymetry = swan_botread(self.mesh_folder+"/"+mesh_bot)


    def load_timesteps(self):
        """ (tentative solution: eventually can check start dates from the s3 bucket itself.
            self.start_year and self.start_month are initially from the input_master.json,
            so all timesteps are always from beginning of earliest data folder,
            even after startdate.txt is updated
        """
        datenum = date2num(
            datetime(int(self.start_year), int(self.start_month), 1),
            units="hours since 1970-01-01 00:00:00.0",
            calendar="gregorian"
        )
        self.timesteps = [
            num2date(
                datenum + t,
                units="hours since 1970-01-01 00:00:00.0",
                calendar="gregorian"
            )
            for t in range(130000) # tentative
        ]
        if os.path.exists("startdate.txt"):
            with open("startdate.txt", "r") as sd:
                self.start_date = sd.readline()
                self.start_year = self.start_date[:4]
                self.start_month = self.start_date[5:7]
                print(f"read startdate.txt: |{self.start_date}|, year: {self.start_year}, month: {self.start_month}")
        else:
            print("no file called \'startdate.txt\'")

        self.num_timesteps = len(self.timesteps)


    def create_nca_input(self):
        """
        Creates an 'nca' structure for the master input template (json) using the swan mesh data.
        No variable data is loaded yet.
        (most of this should already be in input_master.json for the user to specify)
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

        dimensions = dict(
            npe=3,
            nnode=self.num_nodes,
            ntime=self.num_timesteps,
            nelem=self.num_elements
        )
        master_input["nca"]["dimensions"] = dimensions

        variables = dict(
            lat=dict(type="float64", dimensions=["nnode"], units="degrees_north", standard_name="latitude", long_name="latitude"),
            lon=dict(type="float64", dimensions=["nnode"], units="degrees_east", standard_name="longitude", long_name="longitude"),
            elem=dict(type="int32", dimensions=["nelem", "npe"], units="", standard_name="elements", long_name="elements"),
            time=dict(type="float64", dimensions=["ntime"], units="hours since 1970-01-01 00:00:00.0", calendar="gregorian", standard_name="Datetime", long_name="Datetime"),
            b=dict(type="float32", dimensions=["nnode"], units="m", standard_name="Bathymetry", long_name="Bathymetry, m (CGVD28)")
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



    # =====================================
    def load_mat(self, filepath, mfile):
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
        timestep_keys = []

        #print(f" (load_mat) reading from {filepath}")
        mname = re.match(re_mat, mfile).groups()[0] # remove 'test' when ready <------

        # mname: HS
        # mfile: HS.mat
        # filepath: ../BCSWANv5/2004/01/results/HS.mat

        matfile_dict = {}
        matfile_dictx = {}
        matfile_dicty = {}
        #print(mfile)
        try:
            matfile = loadmat(filepath)
            matfile.pop('__header__')
            matfile.pop('__version__')
            matfile.pop('__globals__')

            # k: 'Hsig_20040101_030000' or 'Windx_20040106_120000', etc
            # v: 'array([-22.5, -18.0, -36.7, ...])' for all 177495 nodes
            for k, v in matfile.items():
                v_sq = np.squeeze(v)
                ts = k[-15:]
                timestep = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                timestep_keys.append(timestep)

                if re.match(re_x, k):
                    matfile_dictx.update({k: v_sq})
                elif re.match(re_y, k):
                    matfile_dicty.update({k: v_sq})
                else:
                    matfile_dict.update({k: v_sq})

            if len(matfile_dict) == 0:
                self.matfile1[mname+"X"] = matfile_dictx
                self.matfile2[mname+"Y"] = matfile_dicty
            else:
                self.matfile1[mname] = matfile_dict

            self.start_date = timestep_keys[0]  # first timestep in the month


        except NotImplementedError as e:
            print(f"{e}!")


    def load_spc(self, filepath, sfile):
        """
            spectra files are stored as .txt instead of .mat

            number of timesteps:
                - equal to month, same as mat files.

                number of FACTOR blocks for that timestep:
                    - equal to number of lat+lon nodes (e.g. 22 pairs)

                    create table for that FACTOR block:
                        - columns are each direction (e.g. 36 across)
                        - rows are each frequency (e.g. 34 down)

            performance issue: dictionaries vs numpy arrays? mix? or use one over the other?
        """
        timestep_keys = []

        #print(f" (load_spc) reading from {filepath}")
        sname = re.match(re_spc, sfile).groups()[0]

        s = open(filepath, "r")

        lonlats = []
        afreqs = []
        dirs = []
        num_lonlat, num_afreq, num_dir = 0., 0., 0.

        spcfile = {}

        metadata = True
        while metadata:
            line = s.readline()
            token = line.split()[0]

            if token == "LONLAT":
                num_lonlat = int(s.readline().split()[0])
                for ll in range(num_lonlat):
                    lonlat = s.readline().split()
                    lonlats.append((lonlat[0], lonlat[1]))
                #print("LONLAT ->", num_lonlat)

            elif token == "AFREQ":
                num_afreq = int(s.readline().split()[0])
                for af in range(num_afreq):
                    afreq = s.readline().split()[0]
                    afreqs.append(afreq)
                #print("AFREQ ->", num_lonlat)

            elif token == "NDIR":
                num_dir = int(s.readline().split()[0])
                for d in range(num_dir):
                    dir = s.readline().split()[0]
                    dirs.append(dir)
                #print("NDIR ->", num_dir)

            elif re.match(re_spcdate, token):
                metadata = False

        self.afreqs = afreqs
        self.ndir = dirs
        # some more info between above and tables below

        data = True
        while data:
            token = line.split()[0]
            if re.match(re_spcdate, token):
                timestep = datetime(int(token[:4]), int(token[4:6]), int(token[6:8]), int(token[9:11]))
                timestep_keys.append(timestep)

                lldict = {}
                for lonlat in lonlats:
                    FACTOR = s.readline().split()[0].strip()

                    if FACTOR == "NODATA":  # no data for this lat-long
                        continue

                    factor = s.readline().split()[0].strip()
                    factor = float(factor)

                    fd = np.array(
                        [np.array([int(d) for d in s.readline().split()]) for afreq in range(num_afreq)]
                    )
                    lldict[lonlat] = (factor, fd)
                spcfile[timestep] = lldict
            else:
                print("something went wrong!")
                data = False

            line = s.readline()  # ready next line
            if not line:
                #print("no line!")
                data = False

        s.close()

        self.start_date = timestep_keys[0]  # first timestep in the month
        self.spcfile = spcfile


    def upload_to_cache(self, filetype):
        """
            A NetCDF2D object is created with the finished template and the contents from the node map
                are written to the NetCDF2D object one timestep at a time.
            if 'grid' is true, that means we need to re-load all static data as well as the time steps.

            upload_to_cache then empties both self.matfile1 and self.matfile2 --> {} for the next upload
        """
        netcdf2d = NetCDF2D(self.master_input)

        # grid ('static' data) if reloading (only after session interruption or initial upload)
        # ---------------------------
        if filetype == "grid":
            #print("uploading grid data to cache...")
            #start = time.time()
            netcdf2d.nca.variables["lat"][:] = self.lat
            netcdf2d.nca.variables["lon"][:] = self.lon
            netcdf2d.nca.variables["elem"][:] = self.elements
            netcdf2d.nca.variables["b"][:] = self.bathymetry

            for t, ts in enumerate(self.timesteps):
                netcdf2d.nca.variables["time"][t] = date2num(
                    self.timesteps[t],
                    units="hours since 1970-01-01 00:00:00.0",
                    calendar="gregorian"
                )

            #end = time.time()
            #print("`-- grid upload time:", end-start, "\n")
        # ---------------------------

        elif filetype == "mat":
            MAT_val1, MAT_val2 = "", ""
            variable1, variable2 = "", ""
            MAT_val1 = list(self.matfile1.keys())[0]
            if self.matfile2:
                MAT_val2 = list(self.matfile2.keys())[0]

            for kvar, val in self.temp_var_names.items():
                if val["matfile name"] == MAT_val1:
                    variable1 = kvar
                if val["matfile name"] == MAT_val2 and MAT_val2 != "":
                    variable2 = kvar

            m1start = time.time()
            for kv, n in self.matfile1[MAT_val1].items():
                ts = kv[-15:]
                date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                i = self.timesteps.index(date) # getting the index for that timestep
                try:
                    netcdf2d["s", variable1, i] = n
                except ValueError:
                    raise
            #m1end = time.time()
            #print("            `-- m1 upload time:", m1end - m1start)

            if self.matfile2:
                for kv, n in self.matfile2[MAT_val2].items():
                    ts = kv[-15:]
                    date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                    i = self.timesteps.index(date)
                    try:
                        netcdf2d["s", variable2, i] = n
                    except ValueError:
                        raise
                #m2end = time.time()
                #print("            `-- m2 upload time:", m2end - m1end)

            # erase for the next .mat file in the current month folder
            self.matfile1 = {}
            self.matfile2 = {}

        elif filetype == "spc":
            # TODO
            pass


    def upload_files(self):
        """
        - Pseudocode -

        determine starting month folder
        for each year:
            for each month:
                for each .mat (and .spc) file:
                    load .mat file and first timestep into the Node Map (e.g. HS.mat, line_n.spc)
                    use s3-netcdf package to partition the data and store in the cache
                write latest timestep to 'startdate.txt'


        both mat and spc files will be in the folder, but the mat files are 'prioritized'...
            if spc files are read first, will have timesteps but the number of nodes will be completely different and
            not part of the template

        """
        #start_ = time.time()

        # 'startdate.txt' stores the time step from the month where the upload stopped.
        # a new NodeMap sets self.start_year and self.start_month initially. if reloading, upload mesh data first
        sd = open("startdate.txt", "w+")  # should open each month? each year? or does it matter?
        data_folder = self.data_folder
        search_year = self.start_year == ""
        search_month = self.start_month == ""
        reloading = search_year != ""
        if reloading:
            self.upload_to_cache("grid")

        # print("BCSWANv5")
        for year in sorted(os.listdir(data_folder)):
            if year != 'Mesh' and not year.startswith('.') and os.path.isdir(data_folder + "/" + year):  # avoid '.DS_Store' file
                # print("|--", year)
                if not search_year:  # skip through until year is found
                    if year != self.start_year: continue
                    else: search_year = True
                for month in sorted(os.listdir(data_folder + "/" + year)):
                    if not month.startswith('.') and os.path.isdir(data_folder + "/" + year + "/" + month):
                        # print("    |--", month)
                        if not search_month:  # skip through until month is found
                            if month != self.start_month: continue
                            else: search_month = True
                        results = os.path.join(data_folder, year, month, "results")
                        files = [name for name in os.listdir(results)]  # e.g. ['HS.mat'] or ['HS.mat', 'WIND.mat', ... ]
                        for filename in files:

                            if re.match(re_mat, filename) or re.match(re_mat_test, filename):
                                # start = time.time()
                                # print("        |--", filename)
                                mfilepath = os.path.join(data_folder, year, month, "results", filename)
                                self.load_mat(mfilepath, filename)  # updates self.start_date and current mat files

                                # endloadmat = time.time()
                                # print("            |-- load_mat time:", endloadmat - start)
                                self.upload_to_cache("mat")  # write matfile data to NetCDF2D object

                                # endupload = time.time()
                                # print("            `-- upload_to_cache time:", endupload - endloadmat)

                            elif re.match(re_spc, filename):  # .spc
                                # start = time.time()
                                # print("        |--", filename)
                                sfilepath = os.path.join(data_folder, year, month, "results", filename)
                                self.load_spc(sfilepath, filename) # new spc file updates self.start_date

                                # endloadspc = time.time()
                                # print("            |-- load_spc time:", endloadspc - start)
                                self.upload_to_cache("spc")

                                # endupload = time.time()
                                # print("            `-- upload_to_cache time:", endupload - endloadspc)

                        # store start date
                        sd.seek(0)
                        sd.write(str(self.start_date))
                        sd.truncate()

        sd.close()
        # _end = time.time()
        # print("Total time:", _end-start_)


    # ===================== below not needed anymore?
    def sort_coords(self):
        """ should be after getting coordinates. check if tqdm is installed
        """
        lon_i_sort = np.argsort(self.lon[:])
        lat_i_sort = np.argsort(self.lat[:])

        if TQDM:
            print("\nsorting longitudes...")
            for i in tqdm(range(len(lon_i_sort))):
                self.lon_sort.append(lon_i_sort[i])

            print("\nsorting latitudes...")
            for i in tqdm(range(len(lat_i_sort))):
                self.lat_sort.append(lat_i_sort[i])
        else:
            print("\nsorting latitudes and longitudes...")
            for i in range(len(lon_i_sort)):
                self.lon_sort.append(lon_i_sort[i])
            for i in range(len(lat_i_sort)):
                self.lat_sort.append(lat_i_sort[i])

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

        tg = Triangulation(x, y, self.elements)  # will need to change elements if only a specified region
        plt.figure()  # figsize=(10,10))
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

    def get_node_area(self, box=[0, 0, 0, 0]):
        """ Returns a list of nodes within a boxed region: [west, south, east, north]
        """
        lon_box = [l for i, l in enumerate(self.lon_sort) if (box[0] <= self.lon[i] and self.lon[i] <= box[2])]
        lat_box = [l for i, l in enumerate(self.lat_sort) if (box[1] <= self.lat[i] and self.lat[i] <= box[3])]
        box = list(set(lon_box) & set(lat_box))
        return box

    def get_node_area_delta(self, longitude, latitude, delta):
        """ Returns a list of nodes within a region. The point is at a longitude
            and latitude, with each side at +/- delta degrees from that point
        """
        lon_box = [l for i, l in enumerate(self.lon_sort) if
                   (longitude - delta <= self.lon[i] and self.lon[i] <= longitude + delta)]
        lat_box = [l for i, l in enumerate(self.lat_sort) if
                   (latitude - delta <= self.lat[i] and self.lat[i] <= latitude + delta)]
        box = list(set(lon_box) & set(lat_box))
        return box

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
    # =====================================



def main(*args):
    """
        from the command line:
            $ netcdf-swan.py upload
    """
    # get rid of parentheses and script name
    args = args[0]
    args.pop(0)

    i, start_year, start_month, reloading = 0, "", "", False

    if len(args) > i:
        uploading = re.match(re_upload, args[i]) is not None
        downloading = re.match(re_download, args[i]) is not None
        i += 1
    if len(args) > i and re.match(re_year, args[i]):
        reloading = True
        start_year = args[i]
        i += 1
    if len(args) > i and re.match(re_month, args[i]):
        start_month = args[i]
        i += 1

    if uploading:
        nm = NodeMap() #reloading)
        nm.upload_files() #start_year, start_month) # should be automatic? read from cache? or save to a txt file



    elif downloading:
        nm = NodeMap()
        nm.download()


    print("*** finished")


if __name__ == '__main__':
    args_ = sys.argv
    #print(args_)
    main(args_)