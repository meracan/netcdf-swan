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

    2. loads .mat or .spc files from the 'results' folder for each year and month, OR
        loads one .mat or .spc file from all year and month folders if "t" or "st" group

    3. uploads data to cache and cloud using s3netcdf package

    Basic usage:
    create node map first, then upload the .mat files:

    nmg = NodeMap()
    nmg.upload_files("grid")

    nms = NodeMap()
    nms.upload_files("s")

    nmt = NodeMap()
    nmt.upload_files("t")


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
import pandas as pd
import re

import pprint as pp
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
from s3netcdf import NetCDF2D

import logging
LOG_FORMAT = "%(levelname)s %(asctime)s  --| %(message)s"
logging.basicConfig(
        filename="progress.log",
        level=logging.DEBUG,
        format=LOG_FORMAT
        )
logger = logging.getLogger()

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
re_mesh = re.compile(r"^-?[mM](esh)?$")
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
        self.num_snodes = 0
        self.num_elements = 0
        self.num_timesteps = 0
        self.lon, self.lat = None, None
        self.spc_lon, self.spc_lat = None, None
        self.lon_sort, self.lat_sort = [], []
        self.boundary_markers = None
        self.matfile1 = {}
        self.matfile2 = {}
        self.MAT_names = []
        self.spcfile = {}

        self.elements = None
        self.timesteps = []
        self.start_date = ""
        self.start_year = ""
        self.start_month = ""
        self.bathymetry = []

        with open("../jsons/variable_names.json") as vnames:
            self.var_names = json.load(vnames)
        with open("../jsons/read_master.json") as rt:
            self.read_template = json.load(rt)
        with open("../jsons/input_master.json") as im:
            self.master_input = json.load(im)

            self.mesh_folder = self.master_input["file paths"]["mesh folder"]
            self.data_folder = self.master_input["file paths"]["data folder"]
            self.file_checklist = sorted(self.master_input["file paths"]["file checklist"])
            self.start_year = self.master_input["start year"]
            self.start_month = self.master_input["start month"]
            self.chunk_size = self.master_input["chunk size"]
            self.afreqs = self.master_input["absolute frequencies"]
            self.dirs = self.master_input["directions"]
            self.num_afreqs = len(self.afreqs)
            self.num_dirs = len(self.dirs)

        # for "t" group
        self.curr_variable = ""
        self.chunk_index = 0
        self.curr_nodes1 = np.zeros(shape=(self.chunk_size, 130000))
        self.curr_nodes2 = np.zeros(shape=(self.chunk_size, 130000))

        # below could be optional if this class is used more generally/universally
        self.load_mesh()
        self.load_spc_nodes()
        self.load_timesteps()
        self.create_nca_input()  # timesteps can go in here now that all of them are 'loaded'

        #self.upload_to_cache("grid")

        logger.info(f"***** initialized node map *****")




    def load_mesh(self):
        """
            For the swan methods below, "-1" turns MATLAB's 1 indexing into Python's 0 indexing
        """
        def swan_noderead(nodefile):
            with open(nodefile, 'r') as f:
                V = np.loadtxt(f)
            num_nodes = int(V[0, 0])
            node_order = V[1:, 0].astype(int) - 1
            lon, lat, = V[1:, 1], V[1:, 2]
            boundary_marker = V[1:, 3].astype(int)
            return num_nodes, node_order, lon, lat, boundary_marker

        def swan_eleread(elefile):
            with open(elefile, 'r') as f:
                V = np.loadtxt(f, skiprows=1)
            elements = V[:, 1:] - 1
            elements = elements.astype(int)
            return elements

        def swan_botread(botfile):
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

        self.num_nodes, _, self.lon, self.lat, self.boundary_markers = swan_noderead(self.mesh_folder+"/"+mesh_node)
        self.elements = swan_eleread(self.mesh_folder+"/"+mesh_ele)
        self.num_elements = len(self.elements)
        self.bathymetry = swan_botread(self.mesh_folder+"/"+mesh_bot)


    def load_spc_nodes(self):
        """
            scans through ALL spc files first, just to get spectra nodes.
            assumes node coordinates haven't changed in the span of ~15 years,
             so the first year and month 'results' folder is used as the default.
        """
        spc_lons = []
        spc_lats = []
        snodes = []

        def open_spc(filename):
            s = open(filename, "r")
            searching = True
            while searching:
                line = s.readline()
                token = line.split()[0]
                if token == "LONLAT":
                    searching = False
                    num_lonlat = int(s.readline().split()[0])
                    for ll in range(num_lonlat):
                        lonlat = s.readline().split()
                        spc_lons.append(lonlat[0])
                        spc_lats.append(lonlat[1])
                        snodes.append((lonlat[0], lonlat[1]))
            s.close()

        results = os.path.join(self.data_folder, self.start_year, self.start_month, "results")

        files = sorted([name for name in os.listdir(results)])
        for filename in files:
            if re.match(re_spc, filename):
                sfilepath = os.path.join(results, filename)
                open_spc(sfilepath)

        self.spc_lon = spc_lons
        self.spc_lat = spc_lats
        self.snodes = snodes
        self.num_snodes = len(snodes)


    def load_timesteps(self):
        """
            timesteps are always from beginning of earliest data folder even after start_date.txt is updated.
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
            for t in range(130000)  # ~15+ years
        ]

        self.num_timesteps = len(self.timesteps)


    def create_nca_input(self):
        """
            Creates an 'nca' structure for the master input template (json) using the swan mesh data.
            No variable data is loaded yet.
            Most of this should already be in input_master.json for the user to specify
        """

        master_input = self.master_input

        metadata = dict(
            title="?File description?",
            institution="?Specifies where the original data was produced?",
            source="?The method of production of the original data?",
            history="?Provides an audit trail for modifications to the original data?",
            references="?Published or web-based references that describe the data or methods used to produce it?",
            comment="?Miscellaneous information about the data or methods used to produce it?"
        )
        master_input["metadata"] = metadata

        dimensions = dict(
            npe=3,
            nnode=self.num_nodes,
            ntime=self.num_timesteps,
            nelem=self.num_elements,
            nsnode=self.num_snodes,
            nafreq=self.num_afreqs,
            ndir=self.num_dirs

        )
        master_input["nca"]["dimensions"] = dimensions

        # temporal data
        variablesMat = {}
        for k, v in self.var_names.items():
            var = dict(type="float32", units=v["units"], standard_name=v["standard name"], long_name=v["long name"])
            variablesMat.update({k: var})

        # spectral data (just the one variable, not sure what to call it)
        variablesSpec = {
            "spectra": dict(
                type="float32",
                units="m2/Hz/degr",
                standard_name="VaDens",
                long_name="variance densities in m2/Hz/degr",
                exception_value=-0.9900E+02
            )
        }

        # saving space
        _hr_since = "hours since 1970-01-01 00:00:00.0"
        _dt = "Datetime"
        _abs_freq = "absolute frequencies in Hz"
        _spec_dir = "spectral nautical directions"
        _spec_deg = "spectral nautical directions in degr"

        # groups
        groups = dict(
            elem=dict(dimensions=["nelem", "npe"], variables=dict(
                elem=dict( type="i4", units="", standard_name="elements", long_name="elements") )
            ),
            time=dict(dimensions=["ntime"], variables=dict(
                time=dict( type="f8", units=_hr_since, calendar="gregorian", standard_name=_dt,long_name=_dt) )
            ),
            node=dict(dimensions=["nnode"], variables=dict(
                lat=dict( type="f8", units="degrees_north", standard_name="latitude", long_name="latitude"),
                lon=dict( type="f8", units="degrees_east", standard_name="longitude", long_name="longitude"),
                b=dict( type="f", units="m", standard_name="Bathymetry", long_name="Bathymetry, m (CGVD28)") )
            ),
            snode=dict(dimensions=["nsnode"], variables=dict(
                lat=dict( type="f8", units="degrees_north", standard_name="latitude", long_name="latitude"),
                lon=dict( type="f8", units="degrees_east", standard_name="longitude", long_name="longitude") )
            ),
            afreq = dict(dimensions=["nafreq"], variables=dict(
                afreq=dict( type="f8", units="Hz", standard_name="absolute frequency", long_name=_abs_freq) )
            ),
            dir=dict(dimensions=["ndir"], variables=dict(
                dir=dict( type="f8", units="degrees", standard_name=_spec_dir, long_name=_spec_deg) )
            ),
            s=dict(dimensions=["ntime", "nnode"], variables=variablesMat),
            ss=dict(dimensions=["ntime", "nsnode", "nafreq", "ndir"], variables=variablesSpec),
            t=dict(dimensions=["nnode", "ntime"], variables=variablesMat),
            st=dict(dimensions=["nnode", "ntime", "nafreq", "ndir"], variables=variablesSpec)
        )
        master_input["nca"]["groups"] = groups

        self.master_input = master_input

    # =====================================

    def load_mat_partial(self, filepath, mfile):
        """
            Reads a fraction ("chunk") of node information from the .mat file,
            starting from the node at self.chunk_index.
            Data is stored in self.curr_nodes instead of self.matfile.
        """
        mname = re.match(re_mat, mfile).groups()[0]
        timestep_keys = []
        # mname:    "HS"
        # mfile:    "HS.mat"
        # filepath: "../BCSWANv5/2004/01/results/HS.mat"

        matfile_dict = {}
        matfile_dictx = {}
        matfile_dicty = {}

        matfile = loadmat(filepath)
        matfile.pop('__header__')
        matfile.pop('__version__')
        matfile.pop('__globals__')

        for k, v in matfile.items():
            v_sq = np.squeeze(v)
            if re.match(re_x, k):
                matfile_dictx.update({k: v_sq})
            elif re.match(re_y, k):
                matfile_dicty.update({k: v_sq})
            else:
                matfile_dict.update({k: v_sq})

        if len(matfile_dict) == 0:
            mf1pd = pd.DataFrame.from_dict(matfile_dictx)
            mf2pd = pd.DataFrame.from_dict(matfile_dicty)
            transpose1 = np.array(mf1pd.iloc[ self.chunk_index : self.chunk_index + self.chunk_size ])
            transpose2 = np.array(mf2pd.iloc[ self.chunk_index : self.chunk_index + self.chunk_size ])
            try:
                tds = [
                    self.timesteps.index(t) for t in [
                        datetime(int(k[-15:][:4]), int(k[-15:][4:6]), int(k[-15:][6:8]), int(k[-15:][9:11]))
                        for k in mf1pd.columns
                    ]
                ]
                for i in range(self.chunk_index, self.chunk_index + self.chunk_size):
                    self.curr_nodes1[i][tds] = transpose1[i][:]
                    self.curr_nodes2[i][tds] = transpose2[i][:]
            except ValueError:
                pass

            self.MAT_names = [mname+"X", mname+"Y"]

        else:
            mf1pd = pd.DataFrame.from_dict(matfile_dict)
            transpose = np.array(mf1pd.iloc[self.chunk_index : self.chunk_index + self.chunk_size])
            try:
                tds = [
                    self.timesteps.index(t) for t in [
                        datetime(int(k[-15:][:4]), int(k[-15:][4:6]), int(k[-15:][6:8]), int(k[-15:][9:11]))
                        for k in mf1pd.columns
                    ]
                ]
                for i in range(self.chunk_index, self.chunk_index + self.chunk_size):
                    self.curr_nodes1[i][tds] = transpose[i][:]
            except ValueError:
                pass

            self.MAT_names = [mname]


    def load_mat(self, filepath, mfile):
        """
            'filepath' and 'mfile' could be merged into one?
            Loads one mat file into self.matfile1 and self.matfile2 if there are both X and Y coordinates
            (e.g. WIND.mat), otherwise just stores into self.matfile1.

            The header, version and globals are popped for now. (might need later?)

            'np.squeeze()' is needed to get rid of nested arrays like 'array([[-22.5, -18.0, -36.7, ...]])'

            name stripping is to isolate the date:
                key <- 'abcde_20040101_050000'
                key[-15:] -> '20040101_050000'

            the start date is updated with the first time step of the month at the end
        """
        timestep_keys = []
        mname = re.match(re_mat, mfile).groups()[0]

        # mname:    "HS"
        # mfile:    "HS.mat"
        # filepath: "../BCSWANv5/2004/01/results/HS.mat"

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
            logger.error(f"implementation error in load_mat()  XXX\n")
            raise e


    def load_spc(self, filepath, sfile):
        """
            Aggregates all spc files (stored as .txt instead of .mat) for one month.
            AFREQS and NDIRS are already hardcoded as constants.
            Each datum in the block is multiplied with that block's FACTOR.

            number of timesteps:
                - equal to month, same as mat files.

                number of FACTOR blocks for that timestep:
                    - equal to number of lat+lon nodes (e.g. 22 pairs)

                    create table for that FACTOR block:
                        - columns are each direction "dir" (e.g. 36 across)
                        - rows are each frequency "afreq" (e.g. 34 down)
        """
        timestep_keys = []
        sname = re.match(re_spc, sfile).groups()[0]

        s = open(filepath, "r")

        lonlats = []
        num_lonlat, num_afreq, num_dir = 0., len(self.afreqs), len(self.dirs)

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
            elif re.match(re_spcdate, token):
                metadata = False

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
                    fd_block = np.array(
                        [np.array([float(int(d)*factor) for d in s.readline().split()]) for afreq in range(num_afreq)]
                    )
                    lldict[lonlat] = fd_block
                spcfile[timestep] = lldict
            else:
                logger.warning(f"date mismatch in load_spc(). reading from {sfile} stopped    ???")
                data = False
            line = s.readline()  # ready next line
            if not line:
                # print("no line!")
                data = False

        s.close()

        self.start_date = timestep_keys[0]  # save first timestep in the month
        self.spcfile = spcfile


    def upload_to_cache(self, group_type="s"):
        """
            Writing to NetCDF2D:
            A NetCDF2D object is created with the finished template and the contents from the node map
                are written to the NetCDF2D object one timestep at a time if group "s" or "ss",
                or one or more nodes at a time if group "t" or "st".
            The grid ('static' data) only needs to be uploaded once.

            Groups "t" and "st" use pandas to extract the 'transpose' of data table information.

            upload_to_cache then empties (--> {}) the data for the next upload.
        """
        netcdf2d = NetCDF2D(self.master_input)


        if group_type == "grid":
            logger.info("uploading grid data to cache...")

            # dont need?
            timesteps = np.array([
                date2num(self.timesteps[t], units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")
                for t in range(self.num_timesteps)
            ])

            netcdf2d["time", "time"] = self.timesteps
            netcdf2d["node", "lat"] = self.lat
            netcdf2d["node", "lon"] = self.lon
            netcdf2d["node", "b"] = self.bathymetry
            netcdf2d["snode", "lat"] = self.spc_lat
            netcdf2d["snode", "lon"] = self.spc_lon
            netcdf2d["elem", "elem"] = self.elements
            netcdf2d["afreq", "afreq"] = self.afreqs
            netcdf2d["dir", "dir"] = self.dirs


        elif group_type == "s":
            MAT1, MAT2 = "", ""
            variable1, variable2 = "", ""
            MAT1 = list(self.matfile1.keys())[0]
            if self.matfile2:
                MAT2 = list(self.matfile2.keys())[0]
            for kvar, val in self.var_names.items():
                if val["matfile name"] == MAT1:
                    variable1 = kvar
                if val["matfile name"] == MAT2 and MAT2 != "":
                    variable2 = kvar

            for kv, n in self.matfile1[MAT1].items():
                ts = kv[-15:]
                date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                try:
                    t_index = self.timesteps.index(date)
                    netcdf2d["s", variable1, t_index] = n
                except ValueError as ve:
                    logger.warning(f"value error: \"{ve}\" while uploading {MAT1} to cache/netcdf2d object   ???")


            if self.matfile2:
                for kv, n in self.matfile2[MAT2].items():
                    ts = kv[-15:]
                    date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                    t_index = self.timesteps.index(date)
                    try:
                        netcdf2d["s", variable2, t_index] = n
                    except ValueError as ve:
                        logger.warning(f"value error: \"{ve}\" while uploading {MAT2} to cache/netcdf2d object   ???")


            # erase for the next .mat file in the current month folder
            self.matfile1 = {}
            self.matfile2 = {}


        elif group_type == "t":
            if len(self.MAT_names) > 0:
                MAT1, MAT2 = self.MAT_names[0], ""
                if len(self.MAT_names) == 2:
                    MAT2 = self.MAT_names[1]
                variable1, variable2 = "", ""
                for kvar, val in self.var_names.items():
                    if val["matfile name"] == MAT1:
                        variable1 = kvar
                    if val["matfile name"] == MAT2 and MAT2 != "":
                        variable2 = kvar

                node_start, node_end = self.chunk_index, self.chunk_index + self.chunk_size
                netcdf2d["t", variable1, node_start:node_end, :] = self.curr_nodes1[:][:]

                if len(self.MAT_names) == 2:
                    node_start, node_end = self.chunk_index, self.chunk_index + self.chunk_size
                    netcdf2d["t", variable2, node_start:node_end, :] = self.curr_nodes2[:][:]

            # erase for the next chunk of nodes
            self.MAT_names = []
            self.curr_nodes1 = np.zeros(shape=(self.chunk_size, 130000))
            self.curr_nodes2 = np.zeros(shape=(self.chunk_size, 130000))



        elif group_type == "ss":
            for date, snodes in self.spcfile.items():
                t_index = self.timesteps.index(date)
                for snode, afreqs in snodes.items():  # snode will have an index
                    sn = self.snodes.index((snode[0], snode[1]))
                    netcdf2d["ss", "spectra", t_index, sn] = afreqs

            # erase for the next .spc file in the current month folder
            self.spcfile = {}


        elif group_type == "st":
            # TODO
            pass



    def upload_files(self, group="s"):
        if group == "grid" or group == "mesh":
            self.upload_to_cache("grid")
        if group == "s" or group == "ss":
            self.upload_files_s()
        elif group == "t" or group == "st":
            self.upload_files_t()

    def upload_files_t(self):
        """
        Scans through all years and months for 1 of the results files (e.g. HS.mat) from a checklist.
        Pops from the checklist when that variable is finished.

        - Pseudocode -

        get starting variable, if specified
        for each variable ("curr_file") in checklist:
            save current variable in text file in case run is interrupted
            for each chunk of nodes:
                for each year:
                    for each month:
                        for each file in the results folder:
                            if .mat file:
                                load_mat_partial( .mat )
                            else if .spc file:
                                load_spc_partial( .spc )
                if curr_file is .mat:
                    upload_to_cache("t")
                else if curr_file is .spc:
                    upload_to_cache("st")
                ready index for next chunk of nodes
            pop variable from checklist
        """
        start_ = time.time()

        # read file name that was stored in start_file.txt, if after run was interrupted
        with open("../src/start_file.txt", 'a+') as start_file:
            start_file.seek(0)
            current_file = start_file.readline().strip()

        num_chunks = int(self.num_nodes/self.chunk_size)

        last_file = ""
        while len(self.file_checklist) > 0:
            if not current_file:
                current_file = self.file_checklist[0]
            else:
                if current_file != last_file:
                    while current_file != self.file_checklist[0]:
                        last_file = self.file_checklist.pop(0)
                else:
                    current_file = self.file_checklist[0]

            with open("../src/start_file.txt", 'w+') as start_file:
                start_file.seek(0)
                start_file.write(current_file)

            curr_file_missing = False
            search_file_start = time.time()
            for chunk in range(num_chunks):
                chunk_start_time = time.time()
                for year in sorted(os.listdir(self.data_folder)):
                    if year != 'Mesh' and not year.startswith('.') and os.path.isdir(self.data_folder + "/" + year):  # avoid '.DS_Store' file
                        for month in sorted(os.listdir(self.data_folder + "/" + year)):
                            if not month.startswith('.') and os.path.isdir(self.data_folder + "/" + year + "/" + month):
                                results = os.path.join(self.data_folder, year, month, "results")
                                if current_file in [name for name in os.listdir(results)]:
                                    filepath = os.path.join(results, current_file)
                                    if re.match(re_mat, current_file):
                                        self.load_mat_partial(filepath, current_file)
                                    elif re.match(re_spc, current_file):
                                        #self.load_spc_partial(filepath, current_file)  # TODO
                                        pass
                                else:
                                    curr_file_missing = True

                if re.match(re_mat, current_file):
                    self.upload_to_cache("t")
                    ci, cf = self.chunk_index, self.chunk_index+self.chunk_size
                    if not curr_file_missing:
                        time_taken = float("%3.f"%(time.time() - chunk_start_time))
                        logger.info(f" {current_file} nodes {ci}-{cf} uploaded ({time_taken} s)")
                elif re.match(re_spc, current_file):
                    self.upload_to_cache("st")  # TODO

                self.chunk_index += self.chunk_size

            # reset
            self.chunk_index = 0
            last_file = self.file_checklist.pop(0)
            if not curr_file_missing:
                time_taken = float("%3.f"%(time.time() - chunk_start_time))
                logger.info(f" {current_file} finished ({time_taken} s)   ***")
            else:
                logger.warning(f"{current_file} not found    ???")

        total_time = float("%3.f"%(time.time() - start_))
        logger.info(f"********* \"t\" group uploaded ({total_time} s) *********\n")


    def upload_files_s(self):
        """
        - Pseudocode -

        determine starting month folder
        for each year:
            for each month:
                write latest timestep to 'start_date.txt'
                for each .mat (and .spc) file:
                    load .mat file and first timestep into the Node Map (e.g. HS.mat, line_n.spc)
                    use s3-netcdf package to partition the data and store in the cache

        """
        logger.info(f". . . . . . uploading \"s\" group . . . . . .")
        start_ = time.time()

        data_folder = self.data_folder

        # read date that was stored in start_date.txt, if after run was interrupted
        with open("../src/start_date.txt", 'a+') as sd:
            sd.seek(0)
            start_date = sd.readline().strip()
            start_year = start_date[:4]
            start_month = start_date[5:7]

        found_year = (self.start_year == start_year) or (start_year == "")
        found_month = (self.start_month == start_month) or (start_month == "")

        for year in sorted(os.listdir(data_folder)):
            if year != 'Mesh' and not year.startswith('.') and os.path.isdir(data_folder + "/" + year):  # avoid '.DS_Store' file
                if not found_year:  # skip through until year is found
                    if year != str(start_year): continue
                    else: found_year = True
                for month in sorted(os.listdir(data_folder + "/" + year)):
                    if not month.startswith('.') and os.path.isdir(data_folder + "/" + year + "/" + month):
                        if not found_month:  # skip through until month is found
                            if month != str(start_month): continue
                            else: found_month = True

                        start_date = datetime(int(year), int(month), 1)
                        with open("../src/start_date.txt", "w+") as sd:
                            sd.seek(0)
                            sd.write(str(start_date))
                            sd.truncate()

                        results = os.path.join(data_folder, year, month, "results")
                        files = sorted([name for name in os.listdir(results) if name in self.file_checklist])
                        for filename in files:

                            if re.match(re_mat, filename):
                                start_time = time.time()
                                mfilepath = os.path.join(data_folder, year, month, "results", filename)
                                self.load_mat(mfilepath, filename)  # updates self.start_date and current mat files
                                self.upload_to_cache("s")  # write matfile data to NetCDF2D object
                                time_taken = float("%3.f"%(time.time() - start_time))
                                logger.info(f" {year}/{month}/ {filename} uploaded ({time_taken} s)   ***")

                            elif re.match(re_spc, filename):  # .spc should put into separate method
                                start = time.time()
                                sfilepath = os.path.join(data_folder, year, month, "results", filename)
                                self.load_spc(sfilepath, filename) # new spc file updates self.start_date
                                self.upload_to_cache("ss")
                                time_taken = float("%3.f"%(time.time() - start_time))
                                logger.info(f" {year}/{month}/ {filename} uploaded ({time_taken} s)   ***")

        total_time = float("%3.f"%(time.time() - start_))
        logger.info(f"********* \"s\" group uploaded ({total_time} s) *********\n")


    # ===================== below not needed anymore?
    def sort_coords(self):
        """ should be after getting coordinates. check if tqdm is installed
        """
        lon_i_sort = np.argsort(self.lon[:])
        lat_i_sort = np.argsort(self.lat[:])

        if TQDM:
            print("\nsorting longitudes...")
            for i in range(len(lon_i_sort)):
                self.lon_sort.append(lon_i_sort[i])

            print("\nsorting latitudes...")
            for i in range(len(lat_i_sort)):
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

    def ncdump(self, nc_fid, verb=True):
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
            $ python3 netcdf-swan.py upload
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

    #print(uploading, downloading, load_mesh, start_year, start_month)


    if uploading:
        nm = NodeMap()
        nm.upload_files()

    elif downloading:
        nm = NodeMap()
        nm.download()  # needs parameters

    #print("*** finished ***")


if __name__ == '__main__':
    args_ = sys.argv
    #print(args_)
    main(args_)