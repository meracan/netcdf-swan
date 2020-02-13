#!/usr/bin/env python3

"""
    Testing still allows input parameters on the command line so that a year or month can be specified.

    TODO
        change bare except statements (bad habit!)
"""

import sys, os
import json
import datetime
import pprint as pp
from netCDF4 import num2date, date2num
from datetime import datetime
from createnodes import NodeMap
from createnetcdf import create_nca_input, update_timesteps
from s3netcdf import NetCDF2D
from s3netcdf.netcdf2d_func import NetCDFSummary

with open("../filepath_names.json") as fpnames:
    names = json.load(fpnames)

    TEST_SWAN_FOLDER = names['test_data_folder']
    TEST_MESH_FOLDER = names['test_mesh_folder']
    TEST_NCA_READ = names['swanv5']
    NC_FILE_MESH = names['output_nc_mesh']
    NC_FILE = names['output_nc']


with open("../read_master.json") as rm:
    Input_Read = json.load(rm)


with open("../test_variable_names.json") as vnames:
    TEST_TEMPORAL_VARIABLES = json.load(vnames)


def test_main(*args):
    """

    :param args: <folder where swan data is> <year> <month>
    :return:
    """
    # initial test parameters:  $ ../data 2004 01
    data_folder = TEST_SWAN_FOLDER
    mesh = TEST_MESH_FOLDER

    args = args[0] # ?
    print(args)

    get_year, get_month = "", ""

    if len(args) > 1:
        if args[1] != '.':
            data_folder = args[1]
            print(f"(reading from \'{args[1]}\')")
    if len(args) > 2:
        get_year = args[2]
    if len(args) > 3:
        get_month = args[3]

    print(f"\n*** loading mesh folder into node map ***")
    nm_mesh = NodeMap()
    nm_mesh.load_mesh(mesh)
    num_nodes = nm_mesh.num_nodes
    num_elements = nm_mesh.num_elements

    Master_Input = None
    try:
        # print(" Reading summary of nca in cache...")
        # pp.pprint(NetCDFSummary(TEST_NCA_READ))

        print(" Preparing NetCDF2D object...")
        netcdf2d_read = NetCDF2D(Input_Read)    # if available, read mesh/static data so far. needs 'nca' for it to work

        print(" Found nca. Deleting...")
        netcdf2d_read.cache.delete() # comment out to prevent reset

    except:
        print(" Couldn't find/read/delete/create nca.")          # otherwise create a new 'nca' for master

    finally:
        print(" Creating new one...")
        with open("../input_master.json") as im:
            master_input = json.load(im)
        Master_Input = create_nca_input(master_input, nm_mesh)  # nm_mesh.timesteps is empty.

    print(f"*** loading .mat files into node map ***")
    nm_from_mats = NodeMap()

    # from .mat files
    if get_year:
        if get_month:           # .mat files for one month
            results = os.path.join(data_folder, get_year, get_month, "results")
        else:                   # .mat files for all months of one year
            results = os.path.join(data_folder, get_year)
        nm_from_mats.load_mat(results)

    else:                       # .mat files for ALL months from ALL years
        for year in os.listdir(data_folder):
            if year != 'Mesh' and not year.startswith('.'):
                for month in os.listdir(data_folder + "/" + year):
                    results = os.path.join(data_folder, year, month, "results")
                    nm_from_mats.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    print("*** updating timesteps in nca ***")
    Master_Input = update_timesteps(Master_Input, nm_from_mats)

    print("*** initiating NetCDF2D object ***")
    netcdf2d = None
    try:
        netcdf2d = NetCDF2D(Master_Input)
    except:
        raise

    print("*** loading nca data into NetCDF2D... ***")  # separate function ?

    for kvar, val in TEST_TEMPORAL_VARIABLES.items():
        MAT_val = val['matfile name']  # need to detect HS_6hr and WIND_6hr
        if MAT_val in nm_from_mats.matfiles:
            print(f" {MAT_val}")
            for kv, n in nm_from_mats.matfiles[MAT_val].items():
                ts = kv[-15:]
                date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                i = nm_from_mats.timesteps.index(date)
                try:
                    netcdf2d["s", kvar, i] = n
                except ValueError as ve:
                    raise

    print("*** NetCDF2d created and \'shipped\'. Attempting to read it back... ***")
    try:
        print(" Reading summary of nca in cache...")
        try:
            summary = NetCDFSummary(TEST_NCA_READ)
            # pp.pprint(summary)
            print("...worked.")
        except:
            print("couldn't read!")

        print(" Reading NetCDF2D object...\n")
        netcdf2d_read_again = NetCDF2D(Input_Read)
        print("*    significant wave height of nodes in 4th timestep:", netcdf2d_read_again["s", "hs", 3])
        print("*    wind x velocity of nodes in first timestep:", netcdf2d_read_again["s", "u10", 0])
    except:
        raise

    print("\n*** Finished ***")
  

if __name__ == "__main__":
    args = sys.argv
    test_main(args)