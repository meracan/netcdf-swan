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
from createnodes import NodeMap
from createnetcdf import create_nca_input, update_timesteps,\
    print_netcdf
from s3netcdf import NetCDF2D
import numpy as np
from s3netcdf.netcdf2d_func import createNetCDF, NetCDFSummary,\
  createVariables,getChildShape,getMasterShape,parseDescritor,\
  getIndices,getMasterIndices,getPartitions,dataWrapper

with open("../filepath_names.json") as fpnames:
    names = json.load(fpnames)

    TEST_SWAN_FOLDER = names['test_data_folder']
    TEST_MESH_FOLDER = names['test_mesh_folder']
    TEST_NCA_READ = names['swanv5']
    NC_FILE_MESH = names['output_nc_mesh']
    NC_FILE = names['output_nc']


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


def test_create_master_input():  # after deleting. Doesnt create the .nca file, just the .json object

    data_folder = TEST_SWAN_FOLDER
    mesh = TEST_MESH_FOLDER

    nm_mesh = NodeMap()
    nm_mesh.load_mesh(mesh)
    num_nodes = nm_mesh.num_nodes
    num_elements = nm_mesh.num_elements

    with open("../input_master.json") as im:
        master_input = json.load(im)
    Master_Input = create_nca_input(master_input, nm_mesh)  # nm_mesh.timesteps is empty.

    print("created nca input from mesh")
    assert Master_Input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 0, 'nelem': 348364}
    print("dimensions after loading mesh data:", Master_Input['nca']['dimensions'])
    #pp.pprint(Master_Input)


test_year_month = [('2004', '01'), ('2004', None), (None, None)]



@pytest.mark.parametrize("test_year, test_month", test_year_month)
def test_netcdf2d_static_no_temporal_data(test_year, test_month):

    # load static static data into node map and create nca template
    data_folder = TEST_SWAN_FOLDER
    mesh = TEST_MESH_FOLDER
    nm_mesh = NodeMap()
    nm_mesh.load_mesh(mesh)
    num_nodes = nm_mesh.num_nodes
    num_elements = nm_mesh.num_elements
    with open("../input_master.json") as im:
        master_input = json.load(im)

    # nm_mesh.timesteps is empty
    Master_Input = create_nca_input(master_input, nm_mesh)
    pp.pprint(Master_Input["metadata"])
    assert Master_Input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 0, 'nelem': 348364}

    # load temporal data from mat files into node map (could use the same node map for both)
    nm_from_mats = NodeMap()

    data_folder, get_year, get_month = TEST_SWAN_FOLDER, test_year, test_month

    # pytest parameters test below
    if get_year and get_year is not None:
        if get_month and get_month is not None:  # .mat files for one month
            results = os.path.join(data_folder, get_year, get_month, "results")
        else:  # .mat files for all months of one year
            results = os.path.join(data_folder, get_year)
        nm_from_mats.load_mat(results)
    else:  # .mat files for ALL months from ALL years
        for year in os.listdir(data_folder):
            if year != 'Mesh' and not year.startswith('.'):
                for month in os.listdir(data_folder + "/" + year):
                    results = os.path.join(data_folder, year, month, "results")
                    nm_from_mats.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    # update template now that we have timesteps, using the node map
    Master_Input = update_timesteps(Master_Input, nm_from_mats)
    assert Master_Input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 7, 'nelem': 348364}

    # create NetCDF2D object
    netcdf2d_ = NetCDF2D(Master_Input)

    print("NetCDF2D object created: with static data but no temporal data. Metadata?")
    #print(netcdf2d_["metadata"])
    pp.pprint(netcdf2d_.info())


@pytest.mark.parametrize("test_year, test_month, test_var", test_year_month)
def test_netcdf2d_temporal_data(test_year, test_month, var):
    # attempt to get whats in the cache AFTER static data is in
    print("attempting to read NetCDF2D object from cache (has static data but no temporal data)")
    with open("../read_master.json") as rm:
        read_input = json.load(rm)
    netcdf2d = NetCDF2D(read_input)

    # get data from mat files
    nm_from_mats = NodeMap()

    data_folder, get_year, get_month = TEST_SWAN_FOLDER, test_year, test_month
    if get_year and get_year is not None:
        if get_month and get_month is not None:  # .mat files for one month
            results = os.path.join(data_folder, get_year, get_month, "results")
        else:  # .mat files for all months of one year
            results = os.path.join(data_folder, get_year)
        nm_from_mats.load_mat(results)
    else:  # .mat files for ALL months from ALL years
        for year in os.listdir(data_folder):
            if year != 'Mesh' and not year.startswith('.'):
                for month in os.listdir(data_folder + "/" + year):
                    results = os.path.join(data_folder, year, month, "results")
                    nm_from_mats.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    with open("../test_variable_names.json") as vnames:
        test_variable_names = json.load(vnames)

    # load nca data into NetCDF2D
    for kvar, val in test_variable_names.items():
        MAT_val = val['matfile name']  # need to detect HS_6hr and WIND_6hr
        if MAT_val in nm_from_mats.matfiles:
            for kv, n in nm_from_mats.matfiles[MAT_val].items():
                ts = kv[-15:]
                date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                i = nm_from_mats.timesteps.index(date)
                netcdf2d["s", kvar, i] = n


def test_netcdf2d_HS_then_WIND():
    #netcdf2d_temporal_HS()
    #print("HS loaded")
    netcdf2d_temporal_WIND()
    print("WIND loaded")


# do below two tests one after the other
def netcdf2d_temporal_HS():
    # use this test only after nc files have been created.
    # one variable uploaded at a time

    with open("../read_master.json") as rm:
        read_input = json.load(rm)
    netcdf2d = NetCDF2D(read_input)

    # get data from mat files
    results = os.path.join(TEST_SWAN_FOLDER, '2004', '01', "results")
    nm_from_mats = NodeMap()
    nm_from_mats.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    for kv, n in nm_from_mats.matfiles['HS_6hr'].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm_from_mats.timesteps.index(date)
        netcdf2d["s", "hs", i] = n

def netcdf2d_temporal_WIND():
    # use this test only after nc files have been created.
    # one variable uploaded at a time

    with open("../read_master.json") as rm:
        read_input = json.load(rm)
    netcdf2d = NetCDF2D(read_input)

    # get data from mat files
    results = os.path.join(TEST_SWAN_FOLDER, '2004', '01', "results")
    nm_from_mats = NodeMap()
    nm_from_mats.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    for kv, n in nm_from_mats.matfiles['WIND_6hrX'].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm_from_mats.timesteps.index(date)
        netcdf2d["s", "u10", i] = n



@pytest.mark.parametrize("test_year, test_month", test_year_month)
def test_netcdf2d(test_year, test_month):

    # load static static data into node map and create nca template
    data_folder = TEST_SWAN_FOLDER
    mesh = TEST_MESH_FOLDER
    nm_mesh = NodeMap()
    nm_mesh.load_mesh(mesh)
    num_nodes = nm_mesh.num_nodes
    num_elements = nm_mesh.num_elements
    with open("../input_master.json") as im:
        master_input = json.load(im)
    Master_Input = create_nca_input(master_input, nm_mesh)  # nm_mesh.timesteps is empty.
    assert Master_Input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 0, 'nelem': 348364}

    # load temporal data from mat files into node map (could use the same node map for both)
    nm_from_mats = NodeMap()

    data_folder, get_year, get_month = TEST_SWAN_FOLDER, test_year, test_month

    if get_year and get_year is not None:
        if get_month and get_month is not None:  # .mat files for one month
            results = os.path.join(data_folder, get_year, get_month, "results")
        else:  # .mat files for all months of one year
            results = os.path.join(data_folder, get_year)
        nm_from_mats.load_mat(results)
    else:  # .mat files for ALL months from ALL years
        for year in os.listdir(data_folder):
            if year != 'Mesh' and not year.startswith('.'):
                for month in os.listdir(data_folder + "/" + year):
                    results = os.path.join(data_folder, year, month, "results")
                    nm_from_mats.load_mat(results)  # loads everything into node_map.matfiles and node_map.timesteps

    # update template now that we have timesteps, using the node map
    Master_Input = update_timesteps(Master_Input, nm_from_mats)
    print("dimensions: ", Master_Input['nca']['dimensions'])
    #assert Master_Input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 7, 'nelem': 348364}

    netcdf2d = NetCDF2D(Master_Input)

    with open("../test_variable_names.json") as vnames:
        test_variable_names = json.load(vnames)

    for kvar, val in test_variable_names.items():
        MAT_val = val['matfile name']  # need to detect HS_6hr and WIND_6hr
        if MAT_val in nm_from_mats.matfiles:
            for kv, n in nm_from_mats.matfiles[MAT_val].items():
                ts = kv[-15:]
                date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
                i = nm_from_mats.timesteps.index(date)
                netcdf2d["s", kvar, i] = n


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
            7 timesteps, 177945 nodes
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

    nm = NodeMap()
    nm.load_mat(mtresults)
    t, n = len(nm.timesteps), len(list(nm.matfiles[test_mat].values())[0])
    test_dim_mat = [t, n]
    with Dataset(nc, "r") as src_file:
        tt = src_file.variables[var].get_dims()[0].size  # ntime
        nn = src_file.variables[var].get_dims()[1].size  # nnodes
        test_dim_nc = [tt, nn]

    print(test_dim_mat, test_dim_nc, getChildShape(test_dim_mat), getMasterShape(test_dim_nc))
    #print(getMasterIndices(test_dim_mat), getMasterIndices(test_dim_nc))
    np.testing.assert_array_equal(getChildShape(test_dim_mat), test_dim_nc)
    #np.testing.assert_array_equal(getMasterShape(test_dim_nc), test_dim_mat)

    shape1 = [8, 32768]
    # [1, 1, 8, 32768]
    # [[0], [0,1,2, ..., 32765, 32766, 32767]]
    # [[0], [0]]
    # [[0, 1, 2, 3, 4, 5, 6, 7], [0]]

    shape2a = [7, 177945]
    master2a = getMasterShape(shape2a)                              # [4, 1, 2, 177945]
    indices2a_1 = getIndices((0), shape2a)                          # [[0], [0,1,2, ..., 177942, 177943, 177944]]
    indices2a_2 = getIndices((0, 0), shape2a)                       # [[0], [0]]
    indices2a_3 = getIndices((slice(None, None, None),0), shape2a)  # [[0, 1, 2, 3, 4, 5, 6], [0]]

    print("\n", master2a)
    print("\n", indices2a_1)
    print("\n", indices2a_2)
    print("\n", indices2a_3, "\n")
    print("get master indices shape:", getMasterIndices(indices2a_1, shape2a, master2a).shape)
    print("get master indices shape:", getMasterIndices(indices2a_2, shape2a, master2a).shape)
    print("get master indices shape:", getMasterIndices(indices2a_3, shape2a, master2a).shape)
    print(getPartitions(indices2a_1, shape2a, master2a))
    print(getPartitions(indices2a_2, shape2a, master2a))
    print(getPartitions(indices2a_3, shape2a, master2a))
    #np.testing.assert_array_equal(getMasterIndices(indices2a_1, shape2a, master2a).shape, (32768, 4))
    #np.testing.assert_array_equal(getMasterIndices(indices2a_2, shape2a, master2a).shape, (1, 4))
    #np.testing.assert_array_equal(getMasterIndices(indices2a_3, shape2a, master2a).shape, (8, 4))
    # indices,shape,masterShape

    # x = loadmat('WIND.mat') # or something similar. use createnodes
    #


    print("*** finished ***\n")


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

@pytest.mark.parametrize("year, month, day, ptime, pnode, var", test_time_start_end)
def test_download(year, month, day, ptime, pnode, var):
    # test if s3 download works? 'localOnly' is true at the moment

    print("*** initializing NetCDF2D object")
    with open("../input_master.json") as ri:
        read_input = json.load(ri)
    netcdf2d_read = NetCDF2D(read_input)

    year, month, day, ptime, pnode = int(year), int(month), int(day), int(ptime), int(pnode) # CLI
    name = netcdf2d_read.groups["s"].attributes[var]["standard_name"]
    dim = netcdf2d_read.groups["s"]
    #master = netcdf2d_read.variables["time"]

    #print(master)

    pp.pprint(netcdf2d_read.groups["s"].child)
    #pp.pprint(netcdf2d_read.s3.list())# groups["s"].attributes["master"])


    print(f"{name} of nodes in time-partition {ptime} and node-partition {pnode}:", netcdf2d_read["s", var, [ptime, pnode]])


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

    #test_netcdf2d('2004', '01')
    #test_netcdf2d_read_completed_cache()

    #test_compare_shapes('SWANv5_s_u10_3_0.nc', 'WIND_6hrX.mat')
    #for (nc, mt) in test_nc_against_mat:
    #    test_compare_shapes(nc, mt)

    #print("\n\n- Download Test -")
    #test_download(2004, 1, 1, 0, 0, "hs")

    #for (y, m, d, f, t, v) in test_time_start_end:
    #    test_download(y, m, d, f, t, v)
