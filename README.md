# netcdfswan
netcdf-swan takes SWAN (third-generation wave model) data and uses the s3-netcdf package 
to create and upload netcdf files to an Amazon S3 bucket.

# Installation
## From local folder
```bash
git clone https://github.com/meracan/s3-netcdf.git
pip install -e ./s3-netcdf
git clone https://github.com/meracan/netcdf-swan.git
pip install -e ./netcdf-swan
```

## With conda env and testing
```bash
conda create -n s3netcdf python=3.8
conda activate s3netcdf
git clone https://github.com/meracan/s3-netcdf.git
pip install -e ./s3-netcdf
git clone https://github.com/meracan/netcdf-swan.git
pip install -e ./netcdf-swan
```


# Methodology

For each file in the SWAN dataset, netcdfswan.py will store the file's contents in a NodeMap object, 
then create a NetCDF2D object to automatically partition the NodeMap data into one .nca file and multiple .nc files, 
which are stored in a local cache, beside the SWAN dataset. 
The NetCDF2D object then acts as a 'vehicle' to transport the cache contents to and from the Amazon S3 bucket:

```                                                                     
                                                                ┐       ┌───> s3 Cloud
                                                ┌────>  .nc     |       |
                                                ├────>  .nc     |       |   (NetCDF2D)
                                                ├────>  .nc     |       |    
.mat / .spc  ────>  NodeMap  ────>  NetCDF2D ───┤               ├───────┘ 
                                                |               |
                                                └────>  .nca    |
                                                                ┘
SWAN data                                       local cache
```

For the script to work, a temporally sorted folder structure is assumed for the data, with one 'Mesh' folder plus multiple year folders. 
Each year folder contains 12 month folders. 
Each month folder has a 'results' folder, which contains all of the variable data, 
in the form of Matlab (.mat) and spectra (.spc) files, that is associated with the 
nodes and their coordinate information (the 'mesh' or 'grid'):

```
|
├-- netcdfswan.py
└-- SWAN_DATA
    |
    ├-- 2005
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

The directory path for the SWAN data, as well as the start year and start month, 
should be written in 'input_master.json' ahead of time so the script knows where to look when it runs. 
Ideally, the script should be placed beside the data folder. Other path names are also stored in this json file, 
with the rest of the input as indicated in the s3-netcdf documentation.

### NodeMap

The NodeMap class acts as a sort of 'hub' or temporary storage space to hold all of the static grid and meta data 
from a single .mat or .spc file as nc files are being created, and provides some useful methods for plotting, etc. 
(work in progress...) 
The grid data includes coordinate data such as timesteps, latitudes, longitudes, bathymetry, and other triangular mesh information.

#### ._init\_()

Initializing a NodeMap calls three of its methods:
```
load_mesh()
load_timesteps()
create_nca_input()
```

* ```load_mesh()``` 
extracts the Mesh folder's grid data (.ele, .bot, and .node files) and stores 
it in the NodeMap. No variable data is stored at this point.

* ```load_timesteps()``` 
creates a list of timesteps (```self.timesteps```) for a range of ~20 years,
beginning on the first hour of the start year and month (745 hours per month) indicated in input_master.json. 
The starting date is read from 'startdate.txt' and stored in ```self.start_year``` and ```self.start_month```.
(Note: 
This is a tentative solution for when the run is interrupted in the middle of uploading terabytes of data.
All timesteps are always from beginning of earliest data folder, even after startdate.txt is updated.
Start dates from the s3 bucket itself will be checked in a later version.)


* ```create_nca_input()``` 
creates and fills a template for the 'master' netcdf file (.nca) that will be created with NetCDF2D.
This is stored as ```self.master_input```.

#### .upload_files()
``` 
    def upload_files():
        upload_to_cache("grid")
        determine starting year and month
        for each year:
            for each month:
                for each file (.mat or .spc) in 'results':
                    if .mat file:
                        load_mat()
                        upload_to_cache("mat")
                    if .spc file:
                        load_spc()
                        upload_to_cache("spc")
                write latest timestep to 'startdate.txt'
```

```upload_to_cache("grid")``` will reload all static data and time steps in the case of 
a session interruption (e.g. power outage), starting from the date and time saved in startdate.txt.

Then, for each file per 'results' folder,

* ```load_map()``` or ```load_spc``` 
will store the file's contents into the NodeMap object. 
Because the data comes in varying formats, ```loadmat``` from the Scipy library is used to read Matlab files, 
and a simple script using ```readline()``` is used to parse through the .spc files. 
(Note: 
The .mat and .spc files will be unordered in the folder, but the mat files are 'prioritized'.
If an spc files is read first, will still have the same timesteps but these nodes are a different grid 
from the ones in the template. This will be fixed soon.)
        
* ```upload_to_cache("mat")``` or ```upload_to_cache("spc")``` 
will then create a NetCDF2D object using the finished template (```self.master_input```) 
and write to the NetCDF2D object with the current NodeMap information one timestep (hour) at a time.
The NodeMap information is then cleared for the next file upload.

* The earliest timestep for that month is then saved to startdate.txt.




## Usage

### Basic
create node map first, then upload the .mat files:
```python
from netcdfswan import NodeMap

nm = NodeMap()
nm.upload_files()

```

### From command line

To run the netcdf-swan.py application from the command line:
```
$ python3 netcdfswan.py -U
```



### Commands

```python
nm.upload_files()
```


```python
nm.ncdump()
```



## Testing



### TODO:

commands

input_master: change localOnly to False for actual run

upload(): give stats, how much data was transferred/uploaded, which data/dates, etc.
    'verbose' option
    
download(): need? based on specificity, including variable:
    $ netcdf-swan.py download <year> <month> <variable>
    useful for visual manipulation, triangle mesh grids, .nc creation, etc.

option to convert to non-partitioned .nc file

upload_to_cache("spc")

the 3 init methods could be optional if this class is used more generally/universally


```









```


