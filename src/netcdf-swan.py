#!/usr/bin/env python3

"""
    netcdf-swan.py

    Creates and stores partitioned netcdf files from SWAN data using the s3-netcdf package into a local cache.
    The SWAN data is located on a server, which is where this script will be when it runs.
    The script looks for a Mesh folder (path specified by user or in filepath_names.json) before
        trying to upload anything.
    Entire results folders (one per month) are uploaded at a time (no smaller).

    TODO:
        input_master: change localOnly to False for actual run
        write and read to nca as separate functions? CLI parameters or different .py scripts entirely
            also separate netcdf creation from the s3netcdf manipulation
        may shortcut certain steps (node map is only for convenience)
        upload(): give stats, how much data was transferred/uploaded, which data/dates, etc.
            'verbose' option
        download(): need? based on specificity, including variable:
            $ netcdf-swan.py download <year> <month> <variable>
            useful for visual manipulation, triangle mesh grids, .nc creation, etc.
        option to convert to non-partitioned .nc file
        Since the script will be on the server, should provide the option to upload only chunks
            at a time or re-upload certain parts in case something goes wrong.

    swan data -> node map -> .nc + .nca files (NetCDF2D: write) -> upload (NetCDF2D: S3 client)

"""

import sys, os
import re
import json
import datetime
import pprint as pp
from netCDF4 import num2date, date2num
from datetime import datetime
from createnodes import NodeMap
from createnetcdf import *
from s3netcdf import NetCDF2D
from s3netcdf.netcdf2d_func import NetCDFSummary

# change to normal in src version
with open("../filepath_names.json") as fpnames:
    names = json.load(fpnames)

    SWAN_FOLDER = names['data_folder']  # will have to create an absolute path depending on where the data is
    MESH_FOLDER = names['mesh_folder']
    NCA_READ = names['swanv5']
    # NC_FILE_MESH = names['output_nc_mesh']  # can save nc file
    # NC_FILE = names['output_nc']            # can save nc file

with open("../read_master.json") as rt:
    read_template = json.load(rt)
with open("../input_master.json") as it:
    input_template = json.load(it)
with open("../test_variable_names.json") as vnames:  # !!! remove 'test' when finished
    variable_names = json.load(vnames)


re_upload = re.compile(r"-?[uU](pload)?")
re_download = re.compile(r"-?[dD](ownload)?")
re_edit = re.compile(r"-?[eE](dit)?")  # might not use, but we'll see
re_year = re.compile(r"^20[0-9][0-9]$")
re_month = re.compile(r"^(0?[1-9]|1[0-2])$")
# re_day = re.compile(r"^(0?[1-9]|[12][1-9]|3[01])$")
# re_hour = re.compile(r"^(0?[0-9]|1[0-9]|2[0-3])(00)?$")
# re_lat = re.compile(r"^[-]?lat$")
# re_lon = re.compile(r"^[-]?long?$")


def upload(data_folder, year=0, month=0):
    """
    1. Searches for the Mesh folder in the path 'data_folder' and creates a node map to store the Mesh
        folder's .ele, .bot, and .node data.

    2. Meta data, which includes some from the node map, is put into an input template for the master nc file (.nca)

    3. One 'results' folder is in each month folder, which contains SWAN matlab file data. The contents of the .mat
        files (time coordinate x variable value) are loaded into the node map.

    4. The .nca file template is updated with the timesteps extracted from the .mat data.

    5. A NetCDF2D object is created with the finished template and the contents from the node map
        are written to the NetCDF2D object one timestep at a time.
    """
    print("- Upload -")
    print(f"*** writing mesh data (static) into nca template")
    nm = NodeMap()
    mesh = ""
    if os.path.isdir(data_folder):
        for root, dirs, files in os.walk(data_folder):
            if "Mesh" in dirs:
                mesh = os.path.join(root, "Mesh")
                nm.load_mesh(mesh)

    num_nodes = nm.num_nodes
    num_elements = nm.num_elements
    Master_Input = create_nca_input(input_template, nm)  # nm_mesh.timesteps is empty.

    print(f"*** writing mat data (temporal) into nca template")
    if year:
        if month:       # .mat files for one month
            results = os.path.join(data_folder, year, month, "results")
        else:           # .mat files for all months of one year
            results = os.path.join(data_folder, year)
        nm.load_mat(results)
    else:               # .mat files for ALL months from ALL years
        for year in os.listdir(data_folder):
            for month in os.listdir(data_folder + "/" + year):
                results = os.path.join(data_folder, year, month, "results")
                nm.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    print("*** updating time steps in nca template")
    Master_Input = update_timesteps(Master_Input, nm)

    print("*** initializing NetCDF2D object")
    netcdf2d = NetCDF2D(Master_Input)

    print("*** loading nca data into NetCDF2D object")
    for kvar, val in TEMPORAL_VARIABLES.items():
        MAT_val = val['matfile name']
        if MAT_val in nm.matfiles:
            for kv, n in nm.matfiles[MAT_val].items():
                ts = kv[-15:]
                date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                i = nm.timesteps.index(date)
                try:
                    netcdf2d["s", kvar, i] = n
                except ValueError:
                    raise

    print("*** NetCDF2d created and shipped")


def download(filepath=None, year=0, month=0, variable=None):
    """
        ********* TESTING FOR NOW *************

        download...will need lambda functions to do this, but can still read data from nc cache files into node map,
        and do visual manipulation
        (which parameters for the user to specify?)

        TEMPORAL_VARIABLES doesn't have lats, longs, bathymetry or time info etc. If user wants to see
            one of these, then we need a separate name dictionary for those if it is passed in as a variable
        filepath is where to store the downloaded data

        For now, only extracting from local .nca/.nc, not getting from s3
    """
    print("- Download -")
    print("*** initializing NetCDF2D object")
    netcdf2d_read = NetCDF2D(read_template)



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

    # initial test parameters for upload:  $ ../data 2004 01

    args = args[0]  # get rid of parentheses

    i = 0
    uploading = re.match(re_upload, args[i]) is not None
    downloading = re.match(re_download, args[i]) is not None
    data_folder = SWAN_FOLDER
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

    """
    # Downloading
    if len(args) > i and re.match(re_day, args[i]):
        get_day = args[i]
        i += 1
    if len(args) > i and re.match(re_hour, args[i]):
        get_hour = args[i]
    """

    if uploading:
        upload(data_folder, get_year, get_month)

    elif downloading:
        variable = None
        MATFILE_NAMES = [v['matfile name'] for v in TEMPORAL_VARIABLES.values()]
        if args[-1] in TEMPORAL_VARIABLES or args[-1] in MATFILE_NAMES:
            variable = args[-1]

        download()
        # download(data_folder, get_year, get_month, get_day, get_hour, variable) ???

    print("\n*** Finished ***")


if __name__ == "__main__":
    args_ = sys.argv
    main(args_[1:])