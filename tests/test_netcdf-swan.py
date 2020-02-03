#!/usr/bin/env python3

"""
    swan data -> node map -> netcdf -> s3 object -> deploy
    may shortcut some steps above (swan data -> netcdf)

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
import re
from createnodes import NodeMap
from createnetcdf import create_netcdf, write_netcdf_static, write_netcdf_temporal, print_netcdf

TEST_FILE_PATH = "../output/test.nc"
MESH_FOLDER = "../data/Mesh"
XYZ = "../data/XYZ.mat"
TEMP_2004_01_RESULTS = "../data/2004/01/results"
TEMP_2004_01_RESULTS_WIND = "../data/2004/01/results/WIND.mat"

# re_year = re.compile(r"20[0-9][0-9]")
# re_month = re.compile(r"[0-9][0-9]")


def test_main(*args):
    xyz, mesh = XYZ, MESH_FOLDER
    args = args[0]
    data_folder, get_year, get_month = ".", "", ""
    if len(args) > 1:
        data_folder = args[1]
        print(args, args[1])
    if os.path.isdir(data_folder):
        for root, dirs, files in os.walk(data_folder):
            if "XYZ.mat" in files:
                xyz = os.path.join(root, "XYZ.mat")
            if "Mesh" in dirs:
                mesh = os.path.join(root, "Mesh")
    if len(args) > 2:
        get_year = args[2]
    if len(args) > 3:
        get_month = args[3]

    #results = os.path.join(data_folder, year, month, "results")

    #def mat_to_nodes_to_nc(year, month):

    print(f"\n*** Loading node data from Mesh and XYZ ***")
    nm_from_mats = NodeMap()
    nm_from_mats.load_mesh(mesh, xyz)
    nm_from_mats.sort_coords()

    if get_year:
        search = os.path.join(data_folder, get_year)
        if get_month:
            results = os.path.join(data_folder, get_year, get_month, "results")
            # load mat data for that month of that year
            nm_from_mats.load_mat(results)

        else:
            for month in os.listdir(data_folder+"/"+get_year):
                # check if actual month folder? re necessary if so
                results = os.path.join(data_folder, get_year, month, "results")
                # load mat data from ALL months of that year
                nm_from_mats.load_mat(results)
    else:
        print(data_folder, "\n")
        for year in os.listdir(data_folder):
            #print(year)
            for month in os.listdir(data_folder+"/"+year):
                results = os.path.join(data_folder, year, month, "results")
                # load mat data for ALL months from ALL years
                nm_from_mats.load_mat(results)


    print(f"\n*** Creating NetCDF file as \'{TEST_FILE_PATH}\' ***")
    create_netcdf(TEST_FILE_PATH, nm_from_mats)
    write_netcdf_static(TEST_FILE_PATH, nm_from_mats)
    write_netcdf_temporal(TEST_FILE_PATH, nm_from_mats)

    # mat_to_nodes_to_nc(year, month)

    print(f"\n*** Reading created NetCDF file ***")
    nm = NodeMap()
    nm.load_nc(TEST_FILE_PATH)

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
  