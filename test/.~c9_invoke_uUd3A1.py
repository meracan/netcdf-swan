import numpy as np
from datetime import datetime

npe=3
nelem=200
nnode=100
nstation=100
nsnode=20
ntime=14400
nfreq=3
ndir=5

elem=np.arange(nelem*npe,dtype="i4").reshape((nelem,npe))
time=np.datetime64(datetime(2000,1,1))+np.arange((ntime))*np.timedelta64(1, 'h')
lat=np.arange((nnode),dtype="f8")
lon=np.arange((nnode),dtype="f8")
bed=np.arange((nnode),dtype="f4")
slat=np.arange((nstation),dtype="f8")
slon=np.arange((nstation),dtype="f8")
freq=np.arange((nfreq),dtype="f8")
dir=np.arange((ndir),dtype="f8")

nshape=ntime*nnode
shape=(ntime,nnode)
sgroup={
  "u10":np.arange(nshape,dtype="f4").reshape(shape),
  "v10":np.arange(nshape,dtype="f4").reshape(shape),
  "hs":np.arange(nshape,dtype="f4").reshape(shape),
  "tps":np.arange(nshape,dtype="f4").reshape(shape),
  "tmm10":np.arange(nshape,dtype="f4").reshape(shape),
  "tm01":np.arange(nshape,dtype="f4").reshape(shape),
  "tm02":np.arange(nshape,dtype="f4").reshape(shape),
  "pdir":np.arange(nshape,dtype="f4").reshape(shape),
  "dir":np.arange(nshape,dtype="f4").reshape(shape),
  "dspr":np.arange(nshape,dtype="f4").reshape(shape),
  "qp":np.arange(nshape,dtype="f4").reshape(shape),
  "transpx":np.arange(nshape,dtype="f4").reshape(shape),
  "transpy":np.arange(nshape,dtype="f4").reshape(shape),
  }
nshape=ntime*nnode
shape=(nnode,ntime)
tgroup={
  "u10":np.arange(nshape,dtype="f4").reshape(shape),
  "v10":np.arange(nshape,dtype="f4").reshape(shape),
  "hs":np.arange(nshape,dtype="f4").reshape(shape),
  "tps":np.arange(nshape,dtype="f4").reshape(shape),
  "tmm10":np.arange(nshape,dtype="f4").reshape(shape),
  "tm01":np.arange(nshape,dtype="f4").reshape(shape),
  "tm02":np.arange(nshape,dtype="f4").reshape(shape),
  "pdir":np.arange(nshape,dtype="f4").reshape(shape),
  "dir":np.arange(nshape,dtype="f4").reshape(shape),
  "dspr":np.arange(nshape,dtype="f4").reshape(shape),
  "qp":np.arange(nshape,dtype="f4").reshape(shape),
  "transpx":np.arange(nshape,dtype="f4").reshape(shape),
  "transpy":np.arange(nshape,dtype="f4").reshape(shape),
  }

nshape=nstation*nsnode*ntime*nfreq*ndir
shape=(nstation,nsnode,ntime,nfreq,ndir)
spcgroup={
  "spectra":np.arange(nshape,dtype="f8").reshape(shape)
}


# Nee








