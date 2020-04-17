#!/usr/bin/env python3

"""
    TODO
        download/reading of cache:
            specify timesteps
            specify latitude/longitude area
            comparing shapes of matfile arrays vs nc file arrays
"""

import sys, os, time
import copy
import json
import re
import datetime
import inspect
import pytest
import pprint as pp
from scipy.io import loadmat, savemat
import pandas as pd
from tqdm import tqdm
from netCDF4 import Dataset
from netCDF4 import num2date, date2num
from datetime import datetime
from netcdfswan import NodeMap
from tests import test_spectra
from s3netcdf import NetCDF2D
import numpy as np
from s3netcdf.netcdf2d_func import createNetCDF, NetCDFSummary,\
  createVariables,getChildShape,getMasterShape,parseDescritor,\
  getIndices,getMasterIndices,getPartitions,dataWrapper

re_mat_test = re.compile(r"(WIND_6hr|HS_6hr|HS_)[.]mat")  # separate test from actual?
re_mat = re.compile(r"(WIND|HS|TPS|TMM10|TM01|TM02|PDIR|DIR|DSPR|QP|TRANSP)[.]mat")
re_x = re.compile(r"^.*_x_")
re_y = re.compile(r"^.*_y_")

#  Incremental testing, each one adds a step

def ttest_read_from_cache():
    with open("../read_master.json") as rm:
        input_read = json.load(rm)
    try:
        netcdf2d_read = NetCDF2D(input_read)
    except:
        print("couldn't read cache")


def test_delete_cache():
    # only works if data has been put in.
    # can't delete after cache has already been deleted
    with open("../jsons/read_master.json") as rm:
        input_read = json.load(rm)
    try:
        netcdf2d_read = NetCDF2D(input_read)
        netcdf2d_read.cache.delete()
        print("testing: cache deleted")
    except:
        print("testing: couldn't delete cache")

    try:
        os.remove("../src/start_date.txt")
        print("testing: start_date.txt deleted")
    except:
        print("testing: cannot find start_date.txt")

    try:
        os.remove("../src/start_file.txt")
        print("testing: start_file.txt deleted")
    except:
        print("testing: cannot find start_file.txt")



test_year_month = [('2004', '01'), ('2004', None), (None, None)]


@pytest.mark.parametrize("test_year, test_month", test_year_month)
def ttest_netcdf2d_static_no_temporal_data(test_year, test_month):

    # load static static data into node map and create nca template
    nm_mesh = NodeMap()

    assert nm_mesh.master_input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 0, 'nelem': 348364}

    # load temporal data from mat files into node map
    # could use the same node map for both, but we can try to see how it might be split up in two parts.
    # maybe thats where the metadata went?

    nm_from_mats = NodeMap()
    nm_from_mats.load_mats(test_year, test_month)
    print(nm_from_mats.master_input['nca']['dimensions'])
    #assert nm_from_mats.master_input['nca']['dimensions'] == {'npe': 3, 'nnode': 177945, 'ntime': 6, 'nelem': 348364}


    netcdf2d_ = NetCDF2D(nm_from_mats.master_input)

    print("NetCDF2D object created with static data but no temporal data. Where did the metadata go?")
    # pp.pprint(netcdf2d_.info())


def ttest_netcdf2d_HS_then_WIND():

    netcdf2d_temporal_HS()
    print("HS loaded")
    #netcdf2d_temporal_WIND()
    #print("WIND loaded")


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
    netcdf2d = NetCDF2D(nm_from_mats.master_input)

    print("*** loading nca data into NetCDF2D object")
    for kv, n in tqdm(nm_from_mats.matfiles['HS_6hr'].items()):
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
    netcdf2d = NetCDF2D(nm_from_mats.master_input)

    print("*** loading nca data into NetCDF2D object")
    for kv, n in nm_from_mats.matfiles['WIND_6hrX'].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm_from_mats.timesteps.index(date)
        netcdf2d["s", "u10", i] = n

    print("*** NetCDF2d created and shipped")
    # ----------------------





def test_upload_to_cache_auto():

    nm = NodeMap()

    # need grid first
    nm.upload_to_cache("grid")

    nm.load_mat("../BCSWANv5/2004/01/results/HS.mat", "HS.mat")
    print("matfile:", nm.matfile1)

    # upload to cache clears the matfile contents.
    nm.upload_to_cache("mat")
    print("matfile emptied:", nm.matfile1)

    # load the next month
    nm.load_mat("../BCSWANv5/2004/02/results/HS.mat", "HS.mat")
    print("matfile:", nm.matfile1)

    nm.upload_to_cache("mat")
    print("matfile emptied:", nm.matfile1)


def test_upload_to_cache_manual():

    nm = NodeMap()
    nm.upload_to_cache("grid")
    nm.load_mat("../BCSWANv5/2004/01/results/no_test/HS.mat", "HS.mat") # manual?



    # overriding upload_to_cache for first mat file:
    # ---------------------------------
    netcdf2d = NetCDF2D(nm.master_input)

    MAT_val1, MAT_val2 = "", ""
    variable1, variable2 = "", ""
    MAT_val1 = list(nm.matfile1.keys())[0]
    if nm.matfile2:
        MAT_val2 = list(nm.matfile2.keys())[0]

    for kvar, val in nm.temp_var_names.items():
        if val["matfile name"] == MAT_val1:
            variable1 = kvar
        if val["matfile name"] == MAT_val2 and MAT_val2 != "":
            variable2 = kvar

    m1start = time.time()
    for kv, n in nm.matfile1[MAT_val1].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm.timesteps.index(date)  # getting the index for that timestep
        try:
            netcdf2d["s", variable1, i] = n

        except ValueError:
            raise
    m1end = time.time()
    print("            `-- m1 upload time:", m1end - m1start)

    if nm.matfile2:
        for kv, n in nm.matfile2[MAT_val2].items():
            ts = kv[-15:]
            date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
            i = nm.timesteps.index(date)
            try:
                netcdf2d["s", variable2, i] = n
            except ValueError:
                raise
        m2end = time.time()
        print("            `-- m2 upload time:", m2end - m1end)

    nm.matfile1 = {}
    nm.matfile2 = {}
    # ---------------------------------


    nm.load_mat("../BCSWANv5/2004/02/results/HS.mat", "HS.mat")

    # overriding upload_to_cache for second matfile:
    # ---------------------------------
    netcdf2d = NetCDF2D(nm.master_input)

    MAT_val1, MAT_val2 = "", ""
    variable1, variable2 = "", ""
    MAT_val1 = list(nm.matfile1.keys())[0]
    if nm.matfile2:
        MAT_val2 = list(nm.matfile2.keys())[0]

    for kvar, val in nm.temp_var_names.items():
        if val["matfile name"] == MAT_val1:
            variable1 = kvar
        if val["matfile name"] == MAT_val2 and MAT_val2 != "":
            variable2 = kvar

    m1start = time.time()
    for kv, n in nm.matfile1[MAT_val1].items():
        ts = kv[-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        i = nm.timesteps.index(date)  # getting the index for that timestep
        try:
            netcdf2d["s", variable1, i] = n
        except ValueError:
            raise
    m1end = time.time()
    print("            `-- m1 upload time:", m1end - m1start)

    if nm.matfile2:
        for kv, n in nm.matfile2[MAT_val2].items():
            ts = kv[-15:]
            date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
            i = nm.timesteps.index(date)
            try:
                netcdf2d["s", variable2, i] = n
            except ValueError:
                raise
        m2end = time.time()
        print("            `-- m2 upload time:", m2end - m1end)

    nm.matfile1 = {}
    nm.matfile2 = {}
    # ---------------------------------


def mat_table(filename):
    """
    creates a pandas dataframe.
    """
    mfile = re.search(re_mat, filename).groups()[0]  # remove 'test' when ready
    #print(filename, "-->", mfile)

    matfile_dictx = {}
    matfile_dicty = {}
    try:
        matfile = loadmat(filename)
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
                matfile_dictx.update({k: v_sq})
    except NotImplementedError as e:
        print(f"{e}!")

    return matfile_dictx, matfile_dicty


def test_upload_to_cache_transp():

    # attempting transpose of data. (Only 1d matfiles, not WIND.mat, etc.)
    nm = NodeMap()

    nm.upload_to_cache("grid")

    m1start = time.time()
    # overriding upload_to_cache for first mat file:
    # ---------------------------------
    nm.matfile1["HS"], _ = mat_table("../BCSWANv5/2004/01/results/notest/HS.mat")
    mf1pd_1 = pd.DataFrame.from_dict(nm.matfile1["HS"])

    netcdf2d = NetCDF2D(nm.master_input)
    timesteps1 = np.array(mf1pd_1.iloc[0:nm.num_nodes])
    timedates = [
        datetime(int(k[-15:][:4]), int(k[-15:][4:6]), int(k[-15:][6:8]), int(k[-15:][9:11]))
        for k in mf1pd_1.columns
    ]
    tds = [nm.timesteps.index(t) for t in timedates]

    netcdf2d["t", "hs", 0:1000, tds] = timesteps1[:1000][:]
    nm.matfile1 = {}
    # ---------------------------------
    m1end = time.time()
    print("time taken:", m1end-m1start, "\n")


    raise
    m2start = time.time()
    # overriding upload_to_cache for second matfile (next month):
    # ---------------------------------
    nm.matfile1["HS"], _ = mat_table("../BCSWANv5/2004/02/results/HS.mat")
    mf1pd_2 = pd.DataFrame.from_dict(nm.matfile1["HS"])

    netcdf2d = NetCDF2D(nm.master_input)

    timesteps2 = np.array(mf1pd_2.iloc[0:nm.num_nodes])
    timedates = [
        datetime(int(k[-15:][:4]), int(k[-15:][4:6]), int(k[-15:][6:8]), int(k[-15:][9:11]))
        for k in mf1pd_2.columns
    ]
    tds = [nm.timesteps.index(t) for t in timedates]

    netcdf2d["t", "hs", 0:10, tds] = timesteps2[:10][:]
    nm.matfile1 = {}
    # ---------------------------------
    m2end = time.time()
    print("time taken:", m2end-m2start, "\n")
    print("total time:", m2end-m1start)





def test_netcdf2d_s():
    nm = NodeMap()
    nm.upload_files("s")

def test_netcdf2d_ss():
    nm = NodeMap()
    nm.upload_files("ss")

def test_netcdf2d_t():
    nm = NodeMap()
    nm.upload_files("t")

def test_netcdf2d_st():
    nm = NodeMap()
    nm.upload_files("st")




"""
@pytest.mark.parametrize("test_year, test_month", test_year_month)
def ttest_netcdf2d_old(test_year, test_month):
    nm = NodeMap()
    nm.load_mats(test_year, test_month)
    nm.upload()
"""


def test_netcdf2d_read_completed_cache():

    print("attempting to read NetCDF2D object from cache")
    with open("../jsons/read_master.json") as rm:
        read_input = json.load(rm)
    netcdf2d = NetCDF2D(read_input)
    #print("Summary:")
    #pp.pprint(netcdf2d_read_only.info())

    print("\nReading NetCDF2D object...\n")

    #print(f"*    (spatial) hs, 5 0:", netcdf2d["s", "hs", 5, 0])
    #print(f"*    (temporal) hs, 0 5:", netcdf2d["t", "hs", 0, 5])

    #print(f"*    (spatial) hs, 0 5:", netcdf2d["s", "hs", 0, 5])
    #print(f"*    (temporal) hs, 5 0:", netcdf2d["t", "hs", 5, 0])

    for n in range(20):
        rnd_time = np.random.randint(745)
        rnd_node = np.random.randint(10000)
        s = netcdf2d["s", "hs", rnd_time, rnd_node]
        t = netcdf2d["t", "hs", rnd_node, rnd_time]
        assert s == t
        #print(s, "==", t, "?")

    """
    print("\n\nstatic data:")
    print("----time")
    print("size:", len(netcdf2d.nca.variables["time"]))
    for i, tt in enumerate(netcdf2d.nca.variables["time"]):
        print(num2date(tt, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian"))
        if i > 20:
            break
    print("----lat")
    print(netcdf2d.nca.variables["lat"][:20])
    print("----lon")
    print(netcdf2d.nca.variables["lon"][:20])
    print("----b")
    print(netcdf2d.nca.variables["b"][:20])
    print("----elem")
    print(netcdf2d.nca.variables["elem"][:20])

    """




















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
def ttest_compare_shapes(test_nc, test_mat):
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



def test_get_timesteps():
    tstart = 298032 # date2num("2004-01-01 00:00:00", units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")
    tend = tstart + 116064 #date2num("2018-01-01 00:00:00", units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")
    #for t in range(0, 116064):
    #    d = num2date(tstart + t, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")

    #    if d.month == 2:
    #        if d.day >= 29 and d.hour == 0:
    #            print(d)
    #            print(date2num(d, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian"))

        #print(type(d))

    ts = [num2date(298032 + t, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian") for t in range(116064)]

    #print(timesteps[0], timesteps[1], timesteps[-2], timesteps[-1])





def test_transp(matfile):

    # transpose everything
    nm = NodeMap()
    #nm.print_netcdf("../s3/SWANv5/ss/SWANv5_ss_spectra_733_1_0_0.nc")
    #print("***")
    #nm.print_netcdf("../s3/SWANv5/s/SWANv5_s_hs_5_0.nc")
    #print(nm.timesteps) # ..., datetime.datetime(2015, 9, 27, 14, 0), ...
    #print("***")


    matfile = '../BCSWANv5/2004/01/results/' + matfile

    # ----------
    mfile = re.search(re_mat, matfile).groups()[0]  # remove 'test' when ready
    # print(filename, "-->", mfile)
    matfile_dictx = {}
    matfile_dicty = {}
    try:
        matfile = loadmat(matfile)
        matfile.pop('__header__')
        matfile.pop('__version__')
        matfile.pop('__globals__')
        #node = []
        #x = [v for k, v in enumerate(matfile.items())]

        #print(x)

        for k, v in matfile.items():
            v_sq = np.squeeze(v)

            if re.match(re_x, k):
                matfile_dictx.update({k: v_sq})
            elif re.match(re_y, k):
                matfile_dicty.update({k: v_sq})
            else:
                matfile_dictx.update({k: v_sq})

    except NotImplementedError as e:
        print(f"{e}!")

    mf1pd = pd.DataFrame.from_dict(matfile_dictx)
    print(mf1pd.head())
    mf2pd = pd.DataFrame.from_dict(matfile_dicty)


    
    #if len(matfile_dict) == 0:
    #    nm.matfile1[mname + "X"] = matfile_dictx
    #    nm.matfile2[mname + "Y"] = matfile_dicty
    #else:
    #    nm.matfile1[mname] = matfile_dict
    

    #return mf1pd, mf2pd





def test_compare_values(matfile, var1, var2=""):


    #matfile = 'data/2004/01/results/'+matfile
    matfile = '../BCSWANv5/2004/01/results/' + matfile

    with open("../jsons/variable_names.json") as vnames:  # !!! remove 'test' when finished
        temp_var_names = json.load(vnames)

    START = datetime(2004, 1, 1, 0, 0)
    iSTART = date2num(START, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")
    #print("START value:", iSTART)

    mt1, mt2 = mat_table(matfile)
    #print(mt1.head())
    #print(mt2.head())

    # timesteps: first, middle, last, and the two steps between them
    #print(mt1.shape)
    t1 = [0, int(mt1.shape[1]//4), int(mt1.shape[1]//2), int(3*mt1.shape[1]//4), mt1.shape[1]-1]
    t2 = [0, int(mt2.shape[1]//4), int(mt2.shape[1]//2), int(3*mt2.shape[1]//4), mt2.shape[1]-1]

    nodes = [0, 1, 2, 10, 8000, 20000, 177000, 177944, -1]
    #print(mt1.columns[0][:-15], ":")
    for t in t1:
        ts = mt1.columns[t][-15:]
        date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
        iDate = date2num(date, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")
        diff = int(abs(iDate - iSTART))
        i_nc = int(diff//2)
        part = int(diff%2)  # either 1 or 2. this is based on current partition shape. might change for other data
        ncfile = f"SWANv5_s_{var1}_{i_nc}_0.nc"
        ncpath = "../s3/SWANv5/s/"+ncfile

        with Dataset(ncpath, "r") as src_file:
            #print(len(src_file.variables[var1][:]))
            #print(src_file.variables[var1][:][:])
            #print(src_file.variables[var1][0])
            #print(src_file.variables[var1][0][:])
            #print(src_file.variables[var1][0][:][part])
            #print(src_file.variables[var1][:][part][node])
            #print(src_file.variables[var1][:][part][node])
            for n in nodes:
                nc_value = src_file.variables[var1][:][part][n]
                mt_value = mt1.iloc[n][mt1.columns[t]]
                assert nc_value == mt_value
                #print(src_file.variables[var1][:][part][n], mt1.iloc[n][mt1.columns[t]]) #
                #print()
                print(f"{mt1.columns[t]}, node {n} ... \t\tncfile: {nc_value}, matfile: {mt_value}")
        print()

    print("--------")
    if var2:
        for t in t2:
            ts = mt2.columns[t][-15:]
            date = datetime(int(ts[:4]), int(ts[4:6]), int(ts[6:8]), int(ts[9:11]))
            iDate = date2num(date, units="hours since 1970-01-01 00:00:00.0", calendar="gregorian")
            diff = int(abs(iDate - iSTART))
            i_nc = int(diff//2)
            part = int(diff%2)  # either 1 or 2
            ncfile = f"SWANv5_s_{var2}_{i_nc}_0.nc"
            ncpath = "../s3/SWANv5/s/"+ncfile
            with Dataset(ncpath, "r") as src_file:
                for n in nodes:
                    nc_value = src_file.variables[var2][:][part][n]
                    mt_value = mt2.iloc[n][mt2.columns[t]]
                    assert nc_value == mt_value
                    print(f"{mt2.columns[t]}, node {n} ... \t\tncfile: {nc_value}, matfile: {mt_value}")
            print()


    #print("--------")



def create_matfiles():
    #filepath = "../BCSWANv5/2004/01/results/HS_B.mat"
    filepath = None
    matfile = loadmat(filepath)

    num_timesteps = 3+7  # 3 headers + 1 week
    years = [
        "2005", "2006", "2007",
        "2008", "2009", "2010", "2011",
        "2012", "2013", "2014", "2015",
        "2016", "2017", "2018"
    ]

    years = ["2004"]
    months = [
        "01", "02", "03", "04",
        "05", "06", "07", "08",
        "09", "10", "11", "12"
    ]

    for year in years:
        for month in months:
            new_matfile = copy.deepcopy(matfile)
            keys = [k for k in new_matfile.keys()]
            for k in keys:
                if k[:4] == "Hsig":
                    new_matfile[k[:5]+year+month+k[11:]] = new_matfile.pop(k)

            # deleting times, to get accuracy
            del_keys = [k for k in new_matfile.keys()]

            # months with only 30 days
            if month in ["04", "06", "09", "11"]:
                for k in del_keys:
                    if k[11:13] == "31":
                        _ = new_matfile.pop(k)

            # february, including leap years
            if month == "02":
                for k in del_keys:
                    if k[11:12] == "3":
                        _ = new_matfile.pop(k)
                if year not in ["2004", "2008", "2012", "2016"]:
                    for k in del_keys:
                        if k[11:13] == "29":
                            _ = new_matfile.pop(k)

            # deleting times again, to save space
            del_keys = [k for k in new_matfile.keys()]

            for i, k in enumerate(del_keys):
                if i >= num_timesteps:
                    _ = new_matfile.pop(k)

            save_path = "../BCSWANv5/" + year + "/" + month + "/results/HS.mat"
            print(save_path)
            #for k in new_matfile.keys():
            #    print(k)

            # savemat(save_path, new_matfile)
        print()













if __name__ == "__main__":
    args = sys.argv
    #main(args)
    #test_read_from_cache()

    #test_get_timesteps()


    # ---------------------------
    #test_delete_cache()
    #test_upload_to_cache_auto()
    #test_netcdf2d_read_completed_cache()


    # ---------------------------
    #test_delete_cache()
    #test_upload_to_cache_manual()
    #test_netcdf2d_read_completed_cache()


    # ---------------------------
    #test_delete_cache()
    #test_upload_to_cache_transp()
    #test_netcdf2d_read_completed_cache()


    # ---------------------------
    #test_delete_cache()
    #test_netcdf2d_s()
    #test_netcdf2d_read_completed_cache()


    # ---------------------------
    #test_delete_cache()
    #test_netcdf2d_t()
    #test_netcdf2d_read_completed_cache()






    #create_matfiles()

    #test_transp("HS.mat")

    #print("======================== did it work?")
    #ttest_netcdf2d_read_completed_cache()

    #ttest_netcdf2d_static_no_temporal_data('2004', '01') # obsolete

    #test_netcdf2d_temporal_data('2004', '01') # obsolete

    #ttest_netcdf2d_HS_then_WIND() #obsolete

    #test_compare_values('WIND_6hr.mat', "u10", "v10")
    #test_compare_values('HS_6hr.mat', "hs")
    #print()

    #test_compare_values("WIND.mat", "u10", "v10")
    #test_compare_values("HS.mat", "hs")
    #test_compare_values("QP.mat", "qp")

    #test_compare_spc_values("", "spectra")

    #test_shapes()

    #test_compare_shapes('SWANv5_s_u10_3_0.nc', 'WIND_6hrX.mat')
    #for (nc, mt) in test_nc_against_mat:
    #    test_compare_shapes(nc, mt)

    #print("\n\n- Download Test -")
    #test_download(2004, 1, 1, 0, 0, "hs")

    #spectra.scan_spectra("../BCSWANv5/2004/01/results/notest/hotspots.spc")
    #spectra.scan_spectra("../BCSWANv5/2004/01/results/notest/line_n.spc")
    #spectra.scan_spectra("../BCSWANv5/2004/01/results/notest/line_s.spc")
    #spectra.scan_spectra("../BCSWANv5/2004/01/results/notest/line_w.spc")
    #----------
    #test_HS_large()

    print("testing: finished")