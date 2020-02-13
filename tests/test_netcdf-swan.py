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

    :param args: <folder where swan data is> <year> <month> <day> <hour>
    :return:
    """
    # initial test parameters:  $ ../data 2004 01
    data_folder = TEST_SWAN_FOLDER
    mesh = TEST_MESH_FOLDER
    args = args[0]
    get_year, get_month, get_day, get_hour, date = "", "", "", "", ""

    if len(args) > 1:
        if args[1] != '.':
            data_folder = args[1]
            print(f"(reading from \'{args[1]}\')")
    if len(args) > 2:
        get_year = args[2]
    if len(args) > 3:
        get_month = args[3]
    if len(args) > 4:
        get_day = args[4]
        date = datetime(int(get_year), int(get_month), int(args[4]))
        # print("date reference created:", date)
    if len(args) > 5:
        get_hour = args[5]
        date = datetime(int(get_year), int(get_month), int(args[4]), int(args[5]))
        print(date)

    print(f"\n*** Loading node data from Mesh and XYZ ***")
    nm_from_mats = NodeMap()
    nm_from_mats.load_mesh(mesh, xyz)
    nm_from_mats.sort_coords()

    # from .mat files
    if get_year:
        if get_month:   # load mat data for one month
            results = os.path.join(data_folder, get_year, get_month, "results")
        else:           # load mat data for all months of one year
            results = os.path.join(data_folder, get_year)
        nm_from_mats.load_mat(results)

    else:               # load mat data for ALL months from ALL years
        for year in os.listdir(data_folder):
            for month in os.listdir(data_folder + "/" + year):
                results = os.path.join(data_folder, year, month, "results")
                nm_from_mats.load_mat(results)

    print(f"\n*** Creating NetCDF file as \'{NC_TEST_FILE}\' from node map ***")
    create_netcdf(NC_TEST_FILE, nm_from_mats)
    write_netcdf_static(NC_TEST_FILE, nm_from_mats, 'day')
    write_netcdf_temporal(NC_TEST_FILE, nm_from_mats)

    print(f"\n*** Reading created NetCDF file (from \'{NC_TEST_FILE}\') ***")
    nm = NodeMap()
    nm.load_nc(NC_TEST_FILE)

    box = nm.get_node_area_delta(-123.3656, 48.4284, 1.0)
    box2 = nm.get_node_area([-124.3656, 47.4284, -123.3656, 48.4284])

    print(box[:10], "...")
    print("number of nodes around victoria (1 degree box):", len(box2))

    print(type(nm.timesteps[-1]))

    # pyplot image
    # nm.tri_plot()
    # nm.update_plot()

    print("\n*** Finished ***")
  

if __name__ == "__main__":
    args = sys.argv
    test_main(args)