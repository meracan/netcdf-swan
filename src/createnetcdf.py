#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import json
import pprint as pp
from datetime import datetime, timedelta
from createnodes import NodeMap
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from s3netcdf import NetCDF2D
from s3netcdf.netcdf2d_func import createNetCDF, NetCDFSummary

# note: pytest will not show progress bar if used in pycharm!
TQDM = True
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    TQDM = False
  
TEST_FILE_PATH = "nc_examples/test.nc"

with open("../variable_names.json") as vnames:
    TEMPORAL_VARIABLES = json.load(vnames)



def create_nca_input(json_template, node_mesh):
    """
    Creates an 'nca' structure for the master input template (json) using swan mesh data.
    No variable data is loaded yet.

    :param filepath:
    :param node_mesh:
    :return: JSON-like structure ready for NetCDF2D object

    """

    master_input = json_template


    # below is only for part of the Input_Master
    metadata = dict(
        title="File description",
        institution="Specifies where the original data was produced",
        source="The method of production of the original data",
        history="Provides an audit trail for modifications to the original data",
        references="Published or web-based references that describe the data or methods used to produce it",
        comment="Miscellaneous information about the data or methods used to produce it"
    )
    master_input["metadata"] = metadata

    # no timesteps in the static nc file, will be updated as more temporal nc files are created
    dimensions = dict(
        npe=3,
        nnode=node_mesh.num_nodes,
        ntime=0,
        nelem=node_mesh.num_elements
    )
    master_input["nca"]["dimensions"] = dimensions

    variables = dict(
        lat=dict(type="float64", dimensions=["nnode"], units="degrees_north", standard_name="latitude", long_name="latitude"),
        lon=dict(type="float64" ,dimensions=["nnode"] ,units="degrees_east" ,standard_name="longitude" ,long_name="longitude"),
        elem=dict(type="int32", dimensions=["nelem", "npe"], units="", standard_name="elements", long_name="elements"),
        time=dict(type="float64", dimensions=["ntime"], units="hours since 1970-01-01 00:00:00.0", calendar="gregorian",
                  standard_name="Datetime", long_name="Datetime"),
        b=dict(type="float32", dimensions=["nnode"], units="m", standard_name = "Bathymetry", long_name = "Bathymetry, m (CGVD28)")
    )
    master_input["nca"]["variables"] = variables

    # temporal data (check for master)
    variables2 = {}
    for k, v in TEMPORAL_VARIABLES.items():
        var=dict(type="float32", units=v["units"], standard_name=v["standard name"], long_name=v["long name"])
        variables2.update({k: var})

    groups = dict(
        s=dict(dimensions=["ntime", "nnode"], variables=variables2),
    )
    master_input["nca"]["groups"] = groups

    return master_input


def update_timesteps(json_template, node_map):
    json_template["nca"]["dimensions"]["ntime"] = len(node_map.timesteps)
    return json_template


def print_netcdf(filepath):
    with Dataset(filepath, "r") as src_file:
        ncdump(src_file)

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
                print('%s:'.format(ncattr)+repr(nc_fid.variables[key].getncattr(ncattr)))
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


if __name__ == "__main__":
    pass
    #create_netcdf(TEST_FILE_PATH)
    #write_netcdf_static(TEST_FILE_PATH)
    #write_netcdf_temporal(TEST_FILE_PATH)
    #print_netcdf(TEST_FILE_PATH)

