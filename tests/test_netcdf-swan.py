#!/usr/bin/env python3

"""
    TODO
        download/reading of cache:
            specify timesteps
            specify latitude/longitude area
            comparing shapes of matfile arrays vs nc file arrays
"""

import sys, os
import json
import datetime
import inspect
import pytest
import pprint as pp
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from datetime import datetime
from netcdfswan import NodeMap
from s3netcdf import NetCDF2D
import numpy as np
from s3netcdf.netcdf2d_func import createNetCDF, NetCDFSummary,\
  createVariables,getChildShape,getMasterShape,parseDescritor,\
  getIndices,getMasterIndices,getPartitions,dataWrapper


#  Incremental testing, each one adds a step

def test_read_from_cache():
    with open("../read_master.json") as rm:
        input_read = json.load(rm)
    try:
        netcdf2d_read = NetCDF2D(input_read)
    except:
        print("couldn't read cache")


def test_delete_cache():
    # only works if data has been put in.
    # can't delete after cache has already been deleted
    with open("../read_master.json") as rm:
        input_read = json.load(rm)
    try:
        netcdf2d_read = NetCDF2D(input_read)
        netcdf2d_read.cache.delete()
        print("cache deleted")
    except:
        print("couldn't delete cache")


test_year_month = [('2004', '01'), ('2004', None), (None, None)]


@pytest.mark.parametrize("test_year, test_month", test_year_month)
def test_netcdf2d_static_no_temporal_data(test_year, test_month):

    # load static static data into node map and create nca template
    nm_mesh = NodeMap()

    assert nm_mesh.master_input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 0, 'nelem': 348364}

    # load temporal data from mat files into node map
    # could use the same node map for both, but we can try to see how it might be split up in two parts.
    # maybe thats where the metadata went?

    nm_from_mats = NodeMap()
    nm_from_mats.load_mats(test_year, test_month)

    assert nm_from_mats.master_input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 7, 'nelem': 348364}

    netcdf2d_ = NetCDF2D(nm_from_mats.master_input)

    print("NetCDF2D object created with static data but no temporal data. Where did the metadata go?")
    # pp.pprint(netcdf2d_.info())


def ttest_netcdf2d_HS_then_WIND():

    netcdf2d_temporal_HS()
    print("HS loaded")
    netcdf2d_temporal_WIND()
    print("WIND loaded")


# do below two tests one after the other
def netcdf2d_temporal_HS():
    # use this test only after nc files have been created.
    # one variable uploaded at a time

    nm_from_mats = NodeMap()
    nm_from_mats.load_mats('2004', '01')  # ('2004', '01', 'hs') ?

    # 'overriding' nm_from_mats.upload()
    # ----------------------
    print("- Upload -")
    print("*** initializing NetCDF2D object")
    netcdf2d = NetCDF2D(nm_from_mats.read_template)

    print("*** loading nca data into NetCDF2D object")
    for kv, n in nm_from_mats.matfiles['HS_6hr'].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm_from_mats.timesteps.index(date)
        netcdf2d["s", "hs", i] = n

    print("*** NetCDF2d created and shipped")
    # ----------------------


def netcdf2d_temporal_WIND():
    # use this test only after nc files have been created.
    # one variable uploaded at a time

    nm_from_mats = NodeMap()
    nm_from_mats.load_mats('2004', '01')  # ('2004', '01', 'u10') ?

    # 'overriding' nm_from_mats.upload()
    # ----------------------
    print("- Upload -")
    print("*** initializing NetCDF2D object")
    netcdf2d = NetCDF2D(nm_from_mats.read_template)

    print("*** loading nca data into NetCDF2D object")
    for kv, n in nm_from_mats.matfiles['WIND_6hrX'].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm_from_mats.timesteps.index(date)
        netcdf2d["s", "u10", i] = n

    print("*** NetCDF2d created and shipped")
    # ----------------------


@pytest.mark.parametrize("test_year, test_month", test_year_month)
def test_netcdf2d(test_year, test_month):
    nm = NodeMap()
    nm.load_mats(test_year, test_month)  # get_var?
    nm.upload()



def test_netcdf2d_read_completed_cache():

    print("attempting to read NetCDF2D object from cache (has both static and temporal data)")
    with open("../read_master.json") as rm:
        read_input = json.load(rm)
    netcdf2d_read_only = NetCDF2D(read_input)
    print("Summary:")
    pp.pprint(netcdf2d_read_only.info())

    print("\nReading NetCDF2D object...\n")
    try:
        print("*    significant wave height of nodes in 4th timestep:", netcdf2d_read_only["s", "hs", 3])
    except:
        pass
    try:
        print("*    wind x velocity of nodes in first timestep:", netcdf2d_read_only["s", "u10", 0])
    except:
        pass
    #print("*** finished ***")



test_nc_against_mat = [
    ('SWANv5_s_hs_0_0.nc', 'HS_6hr.mat'),
    ('SWANv5_s_hs_1_0.nc', 'HS_6hr.mat'),
    ('SWANv5_s_hs_2_0.nc', 'HS_6hr.mat'),
    ('SWANv5_s_u10_0_0.nc', 'WIND_6hrX.mat'),
    ('SWANv5_s_u10_1_0.nc', 'WIND_6hrX.mat'),
    ('SWANv5_s_u10_2_0.nc', 'WIND_6hrX.mat'),
    ('SWANv5_s_u10_3_0.nc', 'WIND_6hrX.mat'),
    ('SWANv5_s_v10_0_0.nc', 'WIND_6hrY.mat'),
    ('SWANv5_s_v10_1_0.nc', 'WIND_6hrY.mat'),
    ('SWANv5_s_v10_2_0.nc', 'WIND_6hrY.mat')
]


@pytest.mark.parametrize("test_nc, test_mat", test_nc_against_mat)
def test_compare_shapes(test_nc, test_mat):
    """
            6 timesteps, 177945 nodes

            Need to compare arrays obtained from mat files with the arrays created from the .nc files.
    """
    nca = '../s3/SWANv5/SWANv5.nca'
    nc = '../s3/SWANv5/s/'+test_nc
    mtresults = 'data/2004/01/results/'
    test_mat = test_mat[:-4]  # erase '.mat' extension
    var = ""
    with open("../test_variable_names.json") as vnames:
        test_variable_names = json.load(vnames)
        for k,v in test_variable_names.items():
            if test_variable_names[k]["matfile name"] == test_mat:
                var = k
                break

    # from .mat files
    nm = NodeMap()
    nm.load_mat(mtresults) # automatically adds timesteps
    t, n = len(nm.timesteps), len(list(nm.matfiles[test_mat].values())[0])
    test_dim_mat = [t, n]

    # from .nc files
    with Dataset(nc, "r") as src_file:
        tt = src_file.variables[var].get_dims()[0].size  # ntime
        nn = src_file.variables[var].get_dims()[1].size  # nnodes
        test_dim_nc = [tt, nn]

    print("test_dim_mat:", test_dim_mat, " test_dim_nc:", test_dim_nc)
    print("child shape:", getChildShape(test_dim_mat), "master shape:", getMasterShape(test_dim_nc))
    np.testing.assert_array_equal(getChildShape(test_dim_mat), test_dim_nc)

    # print("get master indices:", getMasterIndices(indices, shape, masterShape))
    # print("get master indices shape:", getMasterIndices(indices, shape, masterShape).shape)
    # print("get partitions:", getPartitions(indices2a_1, shape2a, master2a))



def test_shapes():
    shape1 = [8, 32768]
    # [1, 1, 8, 32768]
    # [[0], [0,1,2, ..., 32765, 32766, 32767]]
    # [[0], [0]]
    # [[0, 1, 2, 3, 4, 5, 6, 7], [0]]
    #shape6 = [6, 177945]
    #shape31a = [745, 177945]  # January, March, etc. plus the extra hour at midnight for next month
    #shape31b = [744, 177945]  # January, March, etc. without last extra hour
    #shape30a = [721, 177945]  # April, June, etc. plus the extra hour at midnight for next month
    #shape30b = [720, 177945]  # April, June, etc. without extra midnight hour
    #shape29a = [697, 177945]  # February, leap year, plus March 1st at midnight
    #shape29b = [696, 177945]  # February, leap year, no extra hour
    #shape28a = [673, 177945]  # February, extra hour at midnight on March 1st
    #shape28b = [672, 177945]  # February

    shapes = {
        "shape6": [6, 177945],
        "shape31a": [745, 177945],
        "shape31b": [744, 177945],
        "shape30a": [721, 177945],
        "shape30b": [720, 177945],
        "shape29a": [697, 177945],
        "shape29b": [696, 177945],
        "shape28a": [673, 177945],
        "shape28b": [672, 177945]
    }
    for k, v in shapes.items():
        master = getMasterShape(v)
        child = getChildShape(v)
        print(f"{k}:", v)
        print("child:", child)
        print("master:", master)
        print()


"""
test_time_start_end = [
    (2004, 1, 1, 0, 0, "hs"),
    (2004, 1, 1, 1, 0, "hs"),
    (2004, 1, 1, 2, 0, "hs"),
    (2004, 1, 1, 0, 0, "u10"),
    (2004, 1, 1, 1, 0, "u10"),
    (2004, 1, 1, 2, 0, "u10"),
    (2004, 1, 1, 3, 0, "u10"),
    (2004, 1, 1, 0, 0, "v10"),
    (2004, 1, 1, 1, 0, "v10"),
    (2004, 1, 1, 2, 0, "v10"),
    (2004, 1, 1, 3, 0, "v10")
]
"""
test_time_start_end = [
    (2004, 1, 1, 0, 0, "hs"),
    (2004, 1, 1, 1, 0, "hs"),
    (2004, 1, 1, 2, 0, "hs"),
    (2004, 1, 1, 3, 0, "hs"),
    (2004, 1, 1, 4, 0, "hs"),
    (2004, 1, 1, 5, 0, "hs"),
    (2004, 1, 1, 5, 20000, "hs"),
    (2004, 1, 1, 1, 3000, "hs"),
    (2004, 1, 1, 2, 3, "hs"),
    (2004, 1, 1, 0, 0, "u10"),
    (2004, 1, 1, 1, 0, "u10"),
    (2004, 1, 1, 2, 0, "u10"),
    (2004, 1, 1, 3, 0, "u10"),
    (2004, 1, 1, 4, 0, "u10"),
    (2004, 1, 1, 5, 0, "u10"),
    (2004, 1, 1, 5, 0, "u10"),
    (2004, 1, 1, 3, 30, "u10"),
    (2004, 1, 1, 0, 0, "v10"),
    (2004, 1, 1, 1, 0, "v10"),
    (2004, 1, 1, 2, 0, "v10"),
    (2004, 1, 1, 3, 0, "v10"),
    (2004, 1, 1, 4, 0, "v10"),
    (2004, 1, 1, 5, 0, "v10"),
    (2004, 1, 1, 5, 0, "v10"),
    (2004, 1, 1, 3, 20040, "v10"),
]

@pytest.mark.parametrize("year, month, day, ntime, nnode, var", test_time_start_end)
def ttest_download(year, month, day, ntime, nnode, var):
    # test if s3 download works? 'localOnly' is true at the moment

    print("*** initializing NetCDF2D object")
    with open("../input_master.json") as ri:
        read_input = json.load(ri)
    netcdf2d_read = NetCDF2D(read_input)

    year, month, day, ptime, pnode = int(year), int(month), int(day), int(ntime), int(nnode) # CLI
    name = netcdf2d_read.groups["s"].attributes[var]["standard_name"]
    dim = netcdf2d_read.groups["s"]
    #master = netcdf2d_read.variables["time"]

    #print(master)

    # pp.pprint(netcdf2d_read.groups["s"].child) # array([     2, 177945]
    #pp.pprint(netcdf2d_read.s3.list())# groups["s"].attributes["master"])


    #print(f"{name} of nodes in time-partition {ptime} and node-partition {pnode}:", netcdf2d_read["s", var, 5, 3000])
    print(f"{name} of nodes in time {ntime} and node {nnode}:", netcdf2d_read["s", var, ntime, nnode])


    #print(f"{name} of nodes:", netcdf2d_read["s", var, [0,0]])

    # 'groups': {'s': {..., 'variables': {..., 'master':{'dimensions': ['nmaster']

    # get partition division values (which will determine where a specific timestep is located within the cache)

    #x = netcdf2d_read[]

    # ******* test creating visual map


if __name__ == "__main__":
    args = sys.argv
    #main(args)
    #test_read_from_cache()

    #test_delete_cache()

    #test_netcdf2d_static_no_temporal_data('2004', '01')

    #test_netcdf2d_temporal_data('2004', '01')

    #test_netcdf2d_HS_then_WIND()

    test_netcdf2d('2004', '01')
    #test_netcdf2d_read_completed_cache()

    test_shapes()

    test_compare_shapes('SWANv5_s_u10_3_0.nc', 'WIND_6hrX.mat')
    #for (nc, mt) in test_nc_against_mat:
    #    test_compare_shapes(nc, mt)

    #print("\n\n- Download Test -")
    #test_download(2004, 1, 1, 0, 0, "hs")

    #for (y, m, d, f, t, v) in test_time_start_end:
    #    test_download(y, m, d, f, t, v)
