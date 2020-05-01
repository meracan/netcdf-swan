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
Simple reading usage
```python
from netcdfswan import NetCDFSWAN

input={
  "name":"BCSWANv5",
  "cacheLocation":"../s3",
}

swan=NetCDFSWAN(input)

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
Swan results has stored... 

{variable}
- hs
- tp
- etc

`swan['s',{variable},{timestep}]`

`swan['t',{variable},{nodeId}]`

`swan['spc',{variable},{stationId}]`

**Station Id**



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