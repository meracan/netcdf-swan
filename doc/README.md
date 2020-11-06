# UVIC-BCSWAN API Documentation
## NetCDFSWAN
 
The NetCDFSWAN python package reads and stores SWAN results (third-generation wave model) from/to AWS S3 using the netcdf file format.

## Methodology

NetCDFSWAN is a wrapper for the s3-netcdf package. NetCDFSWAN reads files from a SWAN dataset and uses s3-netcdf to create and store the data into partitioned netcdf (.nc) files into a local cache, then take whatever is in the cache and upload those contents to an Amazon S3 bucket: 

```                                                                     
                                                                ┐       ┌───> s3 Cloud
                                                ┌────>  .nc     |       |
                                                ├────>  .nc     |       |   (NetCDF2D)
                                                ├────>  .nc     |       |    
.mat / .spc  ────> NetCDFSWAN ───>  NetCDF2D ───┤               ├───────┘ 
                                                |               |
                                                └────>  .nca    |
                                                                ┘
SWAN data                                       local cache
```

A single master file (.nca) will contain the variable definitions, metadata and indices, and acts as a reference for the data arrays which are stored in the child files (.nc). These files are stored in a local cache, beside the SWAN dataset. 

NetCDFSWAN interprets the SWAN data as having two types: "static" data, which includes coordinate data such as timesteps, latitudes, longitudes, bathymetry, and other triangular mesh information, and "spatial" data, which consists of results from the SWAN model, stored in MATLAB files (.mat) and spectra files (ASCII/text).

For each .mat or .spc file in the SWAN dataset, NetCDFSWAN will read its contents and hold the data in a temporary NetCDF2D object as it gets partitioned.
The NetCDF2D object will then also act as the 'vehicle' to transport the cache contents to and from the Amazon S3 bucket.

A temporally-sorted folder structure is assumed for the SWAN dataset, with one 'Mesh' folder (coordinate information for the nodes) plus multiple year folders. 
Each year folder contains 12 month folders. 
Each month folder has a 'results' folder, which contains all of the variable data associated with the 
spatial and spectra nodes:

```
|
├-- netcdf-swan (package)
|   ├-- src
|   |   └-- netcdfswan.py
|   |
|   └-- BCSWANv5
|       ├-- BCSWANv5.json ("master input")
|       ├-- BCSWANv5.variables.json
|       :
|
├-- SWAN_DATA
|   |
:   ├-- 2005
    ├-- 2006
    ├-- 2007
    :   ├-- 01
        ├-- 02
        ├-- 03
        :   └-- results
                ├-- HS.mat
                ├-- WIND.mat
                ├-- line_n.spc
                ├-- QP.mat
                :

                :
                ├-- hotspots.spc
                └-- TPS.mat
        :
        ├-- 11
        └-- 12
    :
    ├-- 2016
    ├-- 2017
    |
    └-- Mesh
        ├-- .ele
        ├-- .bot
        └-- .node
```
Ideally, netcdf-swan can be placed beside the data folder ("SWAN_DATA"), 
but the path to the data will be referenced in an 'input' json object when reading/writing to s3 (see below).


## Basic writing usage

Simple writing usage. NetCDFSWAN needs an object with a specific format as described in the following section.
The function `prepareInputJSON` is a helper function to create the object.

```python
from netcdfswan import NetCDFSWAN

swanFolder='tests/data'
jsonFile='tests/json/demo.json'
input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder)

swan=NetCDFSWAN(input)

# Upload lat,lon,elem,datetime,freq,dir
swan.uploadStatic()

# Upload variables in the "s" group
swan.uploadS()

# Upload variables in the "t" group
swan.uploadT()

# Upload variables in the "spc" group
swan.uploadSpc()
```

## Writing usage with logger
```python
import logging

swanFolder='tests/data'
jsonFile='tests/json/demo.json'
input=NetCDFSWAN.prepareInputJSON(jsonFile,swanFolder)

logging.basicConfig(
        filename=os.path.join(swanFolder,input['name'],"progress.log"),
        level=logging.DEBUG,
        format="%(levelname)s %(asctime)s  %(message)s"
        )
logger = logging.getLogger()
swan=NetCDFSWAN(input,logger=logger)

...

```

#### Input Json file

The input for creating a master file contains s3 info, metadata, dimensions, partition group, variables, etc.

Metadata attributes are stored in the metadata object.

Dimensions, groups and variables are stored in the nca object.

The input 'master' file from the BC SWAN dataset is provided as an example:
```json
{
  "name":"SWANv5",
  "cacheLocation":"test/output",
  "bucket":"uvic-bcwave",
  "localOnly":true,
  "nca": {
    "metadata":{
      "title":"BCSWANv5",
      "startDate":"2004-01-01T00:00:00",
      "timeStep(h)":1
    },
    "dimensions" : {"npe":3,"nelem":348364,"nnode":177945,"nstation":27,"nsnode":198,"ntime":130000,"nfreq":34,"ndir":36},
    "groups":{
      "elem":{"dimensions":["nelem","npe"],"variables":{
          "elem":{"type":"i4", "units":"" ,"standard_name":"elements" ,"long_name":"elements"}
        }
      },
      "time":{"dimensions":["ntime"],"variables":{
          "time":{"type":"f8", "units":"hours since 1970-01-01 00:00:00.0","calendar":"gregorian" ,"standard_name":"Datetime" ,"long_name":"Datetime"}
        }
      },
      "nodes":{"dimensions":["nnode"],"variables":{
          "lat":{"type":"f8", "units":"degrees_north" ,"standard_name":"latitude" ,"long_name":"latitude"},
          "lon":{"type":"f8", "units":"degrees_east" ,"standard_name":"longitude" ,"long_name":"longitude"},
          "b":{"type":"f4", "units":"m" ,"standard_name":"Bathymetry" ,"long_name":"Bathymetry, m (CGVD28)"}
        }
      },
      "stations":{"dimensions":["nstation","nsnode"],"variables":{
          "lat":{"type":"f8", "units":"degrees_north" ,"standard_name":"latitude" ,"long_name":"latitude"},
          "lon":{"type":"f8", "units":"degrees_east" ,"standard_name":"longitude" ,"long_name":"longitude"}
        }
      },    
      "freq":{"dimensions":["nfreq"],"variables":{
          "freq":{"type":"f8", "units":"Hz" ,"standard_name":"absolute frequency" ,"long_name":"absolute frequencies in Hz"}
        }
      },
      "dir":{"dimensions":["ndir"],"variables":{
          "dir":{"type":"f8", "units":"degrees" ,"standard_name":"spectral nautical directions" ,"long_name":"spectral nautical directions in degr"}
        }
      },      
      "s":{"dimensions":["ntime","nnode"],"variables":"BCSWANv5.variables.json"},
      "t":{"dimensions":["nnode","ntime"],"variables":"BCSWANv5.variables.json"},
      "spc":{"dimensions":["nstation","nsnode","ntime","nfreq", "ndir"],"variables":"BCSWANv5.specvariables.json"}
    }
  }
}
```
#### Input Variable Json file
A list of variables includes type, units and names associated with that parameter. 
The "matfile name" is taken from the array labels in the data itself.
```json
{
  "u10": {
    "type":"f4",
    "units": "m/s", 
    "standard name": "u Wind velocity", 
    "long name": "u Wind velocity, m/s",
    "matfile name": "Windv_x"
  },
  "v10": {
    "type":"f4",
    "units": "m/s", 
    "standard name": "v Wind velocity", 
    "long name": "v Wind velocity, m/s",
    "matfile name": "Windv_y"
  },
  "hs": {
    "type":"f4",
    "units": "m", 
    "standard name": "significant wave height", 
    "long name": "significant wave height (in s)",
    "matfile name": "Hsig"
  },
  "tps": {
    "type":"f4",
    "units": "s", 
    "standard name": "peak period", 
    "long name": "smoothed peak period (in s)",
    "matfile name": "TPsmoo"
  },
  "tmm10": {
    "type":"f4",
    "units": "s", 
    "standard name": "mean absolute wave period", 
    "long name": "mean absolute wave period (in s)",
    "matfile name": "Tm_10"
  },
  "tm01": {
    "type":"f4",
    "units": "s", 
    "standard name": "mean absolute wave period", 
    "long name": "mean absolute wave period (in s)",
    "matfile name": "Tm01"
  },
  "tm02": {
    "type":"f4",
    "units": "s", 
    "standard name": "mean absolute zero-crossing period", 
    "long name": "mean absolute zero-crossing period (in s)",
    "matfile name": "Tm02"
  },
  "pdir": {
    "type":"f4",
    "units": "degrees", 
    "standard name": "peak wave direction", 
    "long name": "peak wave direction in degrees",
    "matfile name": "Pdir"
  },
  "dir": {
    "type":"f4",
    "units": "", 
    "standard name": "mean wave direction", 
    "long name": "mean wave direction (Cartesian or Nautical convention)",
    "matfile name": "Dir"
  },
  "dspr": {
    "type":"f4",
    "units": "degrees", 
    "standard name": "directional wave spread", 
    "long name": "directional spreading of the waves (in degrees)",
    "matfile name": "Dspr"
  },
  "qp": {
    "type":"f4",
    "units": "", 
    "standard name": "peakedness of wave spectrum",
    "long name": "peakedness of the wave spectrum (dimensionless)",
    "matfile name": "Qp"
  },
  "transpx": {
    "type":"f4",
    "units": "m3/s", 
    "standard name": "x transport of energy",
    "long name": "x transport of energy (in W/m or m3/s)",
    "matfile name": "Transp_x"
  },
  "transpy": {
    "type":"f4",
    "units": "m3/s",
    "standard name": "y transport of energy",
    "long name": "y transport of energy (in W/m or m3/s)",
    "matfile name": "Transp_y"
  }
}
```

#### Input Variable Json file
The spectra data object is similar, with only one variable name "spectra".
```json
{
  "spectra": {
    "type":"f8",
    "units": "m2/Hz/degr",
    "standard name": "VaDens",
    "long name": "variance densities in m2/Hz/degr",
    "exception_value":-0.9900E+02,
    "matfile name": "spectra"
	}
}
```
