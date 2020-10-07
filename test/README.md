# UVIC-BCSWAN API - Testing


The BCSWANv5 folder provides the necessary information for reading from 
the s3 bucket where the BC SWAN dataset is stored. 
Testing allows the user to create an example dataset that imitates the UVic's wave model data. 


## Installation
```bash
conda install pytest
```

## Testing
```bash
python3 test/createData.py
python3 test/test_static.py
python3 test/test_netcdfswan.py
```
Testing is recommended in the order above, with the creation of the "demo" data first.

#### createData.py
This will create ~63.8 MB of data in a folder labelled 'output', in the test folder. 
The output folder will include a 'Mesh' folder containing static grid information, 
as well as one year folder (default '2000'), which has 12 months of .mat and .spc files:
```
├-- netcdf-swan (package)
|   ├-- src
:   ├-- BCSWANv5
    ├-- test
    |   ├-- json
    :   |   ├-- demo.json
        |   ├-- demo.spcvariables.json
        |   └-- demo.variables.json
        ├-- output
        |   ├-- 2000
        :   |   ├-- 1
            |   |   └-- results
            |   |       ├-- ??????.spc
            |   |       ├-- ??????.mat
            |   |       ├-- ??????.mat
            |   |       ├-- ??????.spc
            |   |       :            
            |   ├-- 2
            |   :
            |   ├-- 11
            |   └-- 12  
            | 
            └-- Mesh
```
`check_data()` then makes sure that NetCDFSWAN.load() can read the .spc and .mat data 
with the correct dimensions and shapes as specified in `dataTest.py`

#### test_static.py

This will test reading from the demo dataset. 
`test_static.py` will also test `prepareInputJSON` which prepares an input json object to load into NetCDF2D by 
merging the demo jsons together. 
The input object and metadata can be viewed in the 'output' folder as `demo.spectral.json` and `demo.prepared.json`:
```
    ├-- test
    |   ├-- json
    :   |   ├-- demo.json
        |   ├-- demo.spcvariables.json
        |   └-- demo.variables.json
        ├-- output
        |   ├-- 2000
        :   ├-- Mesh
            ├-- demo.spectral.json
            └-- demo.prepared.json
```



#### test_netcdfswan.py

This will test converting and writing the demo data as .nc files into a local cache.
The location of the cache is next to the netcdf-swan package as a default.
A NetCDFSWAN object with minimal read input is then created to test reading
 .nc files from the cache with the correct dimensions and shapes as specified in `dataTest.py`:

```
input={
    "name":"test1",
    "cacheLocation":"../../s3",
    "localOnly":True
  }
  swan = NetCDFSWAN(input)
```