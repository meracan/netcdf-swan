#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import json
from datetime import datetime, timedelta
from createnodes import NodeMap
from netCDF4 import Dataset
from netCDF4 import num2date, date2num

# note: pytest will not show progress bar if used in pycharm!
TQDM = True
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    TQDM = False
  
TEST_FILE_PATH = "nc_examples/test.nc"

with open("../variable_names.json") as vnames:
    TEMPORAL_VARIABLES = json.load(vnames)


def create_netcdf(filepath, node_map):
    """
    Creates a 'blank' netcdf filewith the necesary attributes but without any data.

    latsort and lonsort are for sorting latitudes and longitudes. The assumption is that it will make
    searching for nodes quicker when specifying a region. Probably won't need later on?

    :param filepath:
    :param node_map:
    :return:
    """
    with Dataset(filepath, "w") as src_file:

        # ------ Metadata -------------
        src_file.title = "File description"
        src_file.institution = "Specifies where the original data was produced"
        src_file.source = "The method of production of the original data"
        src_file.history = "Provides an audit trail for modifications to the original data"
        src_file.references = "Published or web-based references that describe the data or methods used to produce it"
        src_file.comment = "Miscellaneous information about the data or methods used to produce it"

        # ------- Dimensions ----------
        ntime = src_file.createDimension("ntime", None) # dimensionless means dynamic
        nnode = src_file.createDimension("nnode", node_map.num_nodes)
        npe = src_file.createDimension("npe", 3)
        nelem = src_file.createDimension("nelem", len(node_map.elements))

        # ---- Static Variables -------
        elem = src_file.createVariable("elem","i4",("nelem","npe",))
        elem.units = ""
        elem.standard_name = "element"
        elem.long_name = "element"

        lon = src_file.createVariable("lon","f8", ("nnode",))
        lon.units = "degrees_east"
        lon.standard_name = "longitude"
        lon.long_name = "longitude"

        lat = src_file.createVariable("lat", "f8", ("nnode",))
        lat.units = "degrees_north"
        lat.standard_name = "latitude"
        lat.long_name = "latitude"

        lonsort = src_file.createVariable("lonsort","i4", ("nnode",))
        lonsort.units = "node number"
        lonsort.standard_name = "node number"
        lonsort.long_name = "node number of longitude"

        latsort = src_file.createVariable("latsort", "i4", ("nnode",))
        latsort.units = "node number"
        latsort.standard_name = "node number"
        latsort.long_name = "node number of latitude"

        b = src_file.createVariable("b", "f8", ("nnode",))
        b.units = "m"
        b.standard_name = "Bathymetry"
        b.long_name = "Bathymetry, m (CGVD28)"

        time = src_file.createVariable("time", "f8", ("ntime",))
        time.units = "hours since 0001-01-01 00:00:00.0"
        time.calendar = "gregorian"
        time.standard_name = "Datetime"
        time.long_name = "Datetime"

        # ------ Temporal Variables ------
        for k,v in TEMPORAL_VARIABLES.items():
            var = src_file.createVariable(k, "f8", ("ntime", "nnode",))
            var.units = v['units']
            var.standard_name = v['standard name']
            var.long_name = v['long name']

        """
        # ------ Spectral Data -----------
        nnode_c = src_file.createDimension("nnode_c", 2)
        nspectra = src_file.createDimension("nspectra", 50)
        nfreq = src_file.createDimension("nfreq", 33)
        ndir = src_file.createDimension("ndir", 36)
        specnodeindex = src_file.createVariable("specnodeindex", "i4", ("nspectra",))
        spectra = src_file.createVariable("spectra", "f4", ("ntime","nspectra","nfreq","ndir",))
        spectra.units = "TODO"
        spectra.standard_name = "TODO"
        spectra.long_name = "Spectra"
        """


def write_netcdf_static(filepath, node_map, timerange='month' ):
    """
    Trying to limit number of timesteps (temporary)
    extract date elements with node_map.timesteps
    change t[9:11] (hours) if necessary
    do we want to enumerate with i, or use dict keys instead?
      e.g. '20040101_020000': 2004-1-1 2:0:0
    :param filepath:
    :param node_map:
    :param timerange:
    :return:
    """

    with Dataset(filepath, "r+") as src_file:
        # Dimensions
        ntime = len(src_file.dimensions['ntime'])
        nnode = len(src_file.dimensions['nnode'])
        npe = len(src_file.dimensions['npe'])
        nelem = len(src_file.dimensions['nelem'])

        # temporary
        if timerange == 'month': timerange = ntime  # <-----
        elif timerange == 'day': timerange = 24     # <-----
        elif timerange == '3 hours': timerange = 3  # <-----
        elif timerange == '1 hour': timerange = 1   # <-----
        elif timerange == '0': timerange = 0        # <-----

        elem = src_file.variables['elem']
        lon = src_file.variables['lon']
        lat = src_file.variables['lat']
        lonsort = src_file.variables['lonsort']
        latsort = src_file.variables['latsort']
        b = src_file.variables['b']
        time = src_file.variables['time']

        elem[:] = np.array(node_map.elements)
        lon[:] = np.array(node_map.lon)
        lat[:] = np.array(node_map.lat)
        lonsort[:] = np.array(node_map.lon_sort)
        latsort[:] = np.array(node_map.lat_sort)
        b[:] = np.array(node_map.bathymetry)
        print(timerange)
        for i, ts in enumerate(node_map.timesteps):
            if i >= timerange: break
            date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
            time[i] = date2num(date, units=time.units, calendar=time.calendar)


def write_netcdf_temporal(filepath, node_map):
    """
    could be a whole separate netCDF file in the end. need to save space
    store time coord as dictionary? might be faster than list/array
    e.g. '20040101_020000': -56.5
    date2num becomes 17557968.0 for 1 Jan 2004
    :param filepath:
    :param node_map:
    :return:
    """

    with Dataset(filepath, "r+") as src_file:
        time = src_file.variables['time']
        if TQDM:
            for k, v in TEMPORAL_VARIABLES.items():
                var = src_file.variables[k]
                mname = v["matfile name"]
                if mname:
                    try:
                        print(f"\n{mname}:")
                        for tt in tqdm(range(len(time))):
                            var[tt, :] = np.array(node_map.matfiles[mname][node_map.timesteps[tt]])
                    except KeyError as ke:
                        print(f"{mname} file not found.")
        else:
            for k, v in TEMPORAL_VARIABLES.items():
                var = src_file.variables[k]
                mname = v["matfile name"]
                if mname:
                    try:
                        for tt in range(len(time)): # <----- specify actual date/time range
                            var[tt, :] = np.array(node_map.matfiles[mname][node_map.timesteps[tt]])
                        print(f"{mname}")
                    except KeyError as ke:
                        #print(f"({mname} file not found)")
                        pass


        """
        # ------ Spectral Data -----------
        nspectra = len(src_file.dimensions["nspectra"])
        nfreq = len(src_file.dimensions["nfreq"])
        ndir = len(src_file.dimensions["ndir"])
        specnodeindex = src_file.variables['specnodeindex']
        spectra = src_file.variables['spectra']
        specnodeindex[:] = np.arange(0, nspectra, dtype=np.int32)
        """


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
    create_netcdf(TEST_FILE_PATH)
    write_netcdf_static(TEST_FILE_PATH)
    write_netcdf_temporal(TEST_FILE_PATH)
    print_netcdf(TEST_FILE_PATH)

