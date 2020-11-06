# UVIC-BCSWAN API
UVIC-BCSWAN API read and store SWAN results (third-generation wave model) from/to AWS S3. 

## Installation
#### Using pip
```bash
git clone https://github.com/meracan/s3-netcdf.git
pip install ./s3-netcdf
git clone https://github.com/meracan/netcdf-swan.git
pip install ./netcdf-swan
```

#### Using conda env
```bash
conda create -n bcswan python=3.8
conda activate bcswan
git clone https://github.com/meracan/s3-netcdf.git
pip install ./s3-netcdf
git clone https://github.com/meracan/netcdf-swan.git
pip install ./netcdf-swan
```

## Usage


Simple reading usage to download UVic's Swan data from the s3 bucket (needs credentials, see "AWS S3 Credentials" contacts below) 
```python
from netcdfswan import NetCDFSWAN

inputFile = "netcdf-swan/BCSWANv5/BCSWANv5.json"
with open(inputFile, "r") as f: 
  inputJson = json.load(f)

swan = NetCDFSWAN(inputJson)

# Read metadata information
print(swan.info())
```


Simple reading usage, from local cache
```python
from netcdfswan import NetCDFSWAN

input = {
  "name":"BCSWANv5",
  "cacheLocation":"../s3",
}

swan = NetCDFSWAN(input)

# Read lat coordinate for all nodes
print(swan['nodes','lat'])

# Read "hs" for the first time step and all nodes
print(swan['s','hs',0])

# Read "hs" for the first node and all time steps
print(swan['t','hs',0])

# Read "spectra" for the first station and all time steps
print(swan['spc','spectra',0])

```


	
#### Variables

The SWAN data stored on S3 is organized into several groups, where each group has one or more variables:

Static data (**group** - _variables_):
- **elem**  - _elem_ 
- **nodes** - _lat, lon, bed_
- **snodes** - _slat, slon, stationid_
- **stations** - _name_
- **freq** - _freq_
- **dir** - _dir_
- **time** - _time_

Spatial/Temporal data:
- **s**/**t** - _hs, tps, tmm10, tm01, tm02, pdir, dir, dspr, qp, u10, v10, transpx, transpy_

Spectra data:
- **spc** - _spectra_

Values from the SWAN data are accessed via the group name, the name of a variable from that group, and an index or range of indices:

swan['s', {variable}, {timestep}, {nodeId}]

swan['t', {variable}, {nodeId}, {timestep}]

swan['spc', {variable}, {snodeId}, {timestep}, {freq}, {dir}]

examples: 

`swan['s', 'transpy', 20]`  y-component of energy transfer for all timesteps at node 20

`swan['s', 'hs', 0:2, 62083]`  significant wave height for first two timesteps at node 62083

`swan['s', 'hs', :, :4]`  all timesteps for first 4 nodes

`swan['t', 'hs', 62083, 0]`  significant wave height at node 62083 for first timestep

`swan['spc', 'spectra', 260, 5000]`  time index 5000 at last spectra node

`swan['spc', 'spectra', 260, 0, :, 35]`  last spectra node for first timestep, all frequencies but only the last nautical direction (index 35)

`swan['nodes', 'lat']`  all latitudes

`swan['stations', 'name']`  names of all spectra node stations



## Doc and development
[Docs](doc/README.md)

## AWS S3 Credentials
Credentials (for example), `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_DEFAULT_REGION` needs to be save in environment variables. 
For more information, check [link](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html).
Please contact `{name}`

### Data License
[License](LICENSE)

### Code License
[License](LICENSE)
