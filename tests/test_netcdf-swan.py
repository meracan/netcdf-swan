#!/usr/bin/env python3

"""
    swan data -> node map -> netcdf -> s3 object -> deploy
    may shortcut some steps above (swan data -> netcdf)

    'filepath_names' will have to be changed if script is in server, e.g. 'BCSWANv5' instead of '../Data'

    1. Searches for XYZ and mesh files in the path directory specified by the user (no arguments or '.'
      means search current directory) and creates a temporary node map.
      should use that location as a reference.
    Arguments after that are date and time coordinates, to search a specific month of a particular year
      for the 'results' folder, e.g.   $ python3 netcdf-swan.py . 2004 01
      will find the 'results' data folder in January of 2004, in the current directory.

    If a year and month are specified, will search the results folder just in that month.
    Otherwise, will look in the results folder for each month of the year if only that year is specified,
    otherwise ALL results folders will be searched.

    2. Reads all .mat files in the results folder and stores them as node data.
    3. Creates and writes a netcdf file (.nc) with the node data

    4. read node data from nc file (for testing)
    5. plot ?


"""
import sys, os
import json
import datetime
from netCDF4 import num2date, date2num
from createnodes import NodeMap
from createnetcdf import create_netcdf, write_netcdf_static, write_netcdf_temporal, print_netcdf

with open("../filepath_names.json") as fpnames:
    names = json.load(fpnames)

    DATA_FOLDER = names['data']
    NC_TEST_FILE = names['output']
    MESH_FOLDER = names['mesh']
    XYZ = names['xyz']


def test_main(*args):
    xyz, mesh = XYZ, MESH_FOLDER
    args = args[0]
    data_folder, get_year, get_month = DATA_FOLDER, "", ""

    if len(args) > 1:
        if args[1] != '.':
            data_folder = args[1]
            print(f"(reading from \'{args[1]}\')")

    if len(args) > 2:
        get_year = args[2]
    if len(args) > 3:
        get_month = args[3]
    # if len(args) > 4:
    if len(args) > 5:
        date = datetime(int(get_year), int(get_month), int(args[4]), int(args[5]))
        print(date)

    print(f"\n*** Loading node data from Mesh and XYZ ***")
    nm_from_mats = NodeMap()
    nm_from_mats.load_mesh(mesh, xyz)
    nm_from_mats.sort_coords()

    results = ""

    if get_year:
        if get_month:   # load mat data for one month
            results = os.path.join(data_folder, get_year, get_month, "results")
        else:           # load mat data for all months of one year
            results = os.path.join(data_folder, get_year)
        nm_from_mats.load_mat(results)

    else:               # load mat data for ALL months from ALL years
        for year in os.listdir(data_folder):
            # print("Year:", year)
            for month in os.listdir(data_folder+"/"+year):
                # print(" Month:", month)
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
  