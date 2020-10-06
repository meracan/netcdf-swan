import numpy as np
from datetime import datetime

npe=3
nelem=20
nnode=10
nstation=27
nsnode=32
ntime=8760
nfreq=3
ndir=5

elem=np.arange(nelem*npe,dtype="i4").reshape((nelem,npe))
time=np.datetime64(datetime(2000,1,1))+np.arange((ntime))*np.timedelta64(1, 'h')
lat=np.arange((nnode),dtype="f8")
lon=np.arange((nnode),dtype="f8")
nodes=np.column_stack((lon,lat))
bed=np.arange((nnode),dtype="f4")
slat=np.arange((nstation),dtype="f8")
slon=np.arange((nstation),dtype="f8")
freq=np.arange((nfreq),dtype="f8")
dir=np.arange((ndir),dtype="f8")

nshape=ntime*nnode
shape=(ntime,nnode)
variables={
  "WIND":{
    "Windv_x":np.arange(nshape,dtype="f4").reshape(shape),
    "Windv_y":np.arange(nshape,dtype="f4").reshape(shape),
    },
  "HS":{"Hsig":np.arange(nshape,dtype="f4").reshape(shape),},
  "DIR":{ "Dir":np.arange(nshape,dtype="f4").reshape(shape),},
  
  "TPS":{"TPsmoo":np.arange(nshape,dtype="f4").reshape(shape),},
  "TMM10":{"Tm_10":np.arange(nshape,dtype="f4").reshape(shape),},
  "TM01":{"Tm01":np.arange(nshape,dtype="f4").reshape(shape),},
  "TM02":{"Tm02":np.arange(nshape,dtype="f4").reshape(shape),},
 
  "PDIR":{"Pdir":np.arange(nshape,dtype="f4").reshape(shape),},
  "DSPR":{"Dspr":np.arange(nshape,dtype="f4").reshape(shape),},
  "QP":{"Qp":np.arange(nshape,dtype="f4").reshape(shape),},
  "TRANSP":{"Transp_x":np.arange(nshape,dtype="f4").reshape(shape),"Transp_y":np.arange(nshape,dtype="f4").reshape(shape),}
  
  }


nshape=nsnode*ntime*nfreq*ndir
shape=(nsnode,ntime,nfreq,ndir)

spcgroup={
  "spectra":(np.arange(nshape,dtype="f8")).reshape(shape)
}

stations={
    "beverly": 1,
    "brooks": 1,
    "c_dixon": 1,
    "c_eliz": 1,
    "campbell": 1,
    "e_dell": 1,
    "hotspots": 2,
    "line_n": 2,
    "line_w": 3,
    "m_nomad": 1,
    "n_hecat": 1,
    "ne_isle": 1,
    "neah": 2,
    "p_renf": 1,
    "perouse": 1,
    "s_hecat": 1,
    "s_morsb": 1,
    "s_nomad": 1,
    "sombrio": 1,
    "sooke": 1,
    "tarbotn": 1,
    "tillamk": 1,
    "tofino": 1,
    "w_dixon": 1,
    "w_morsb": 1,
    "w_otter": 1,
    "w_washn": 1
}
# Create lat lng for each station
isnode=0
for i,vname in enumerate(stations):
  c=np.array([[1.0,1.0]])
  nsnodes=stations[vname]
  sIndex=isnode
  eIndex=isnode+nsnodes
  latlon=np.zeros((nsnodes,2))
  latlon[:,0]=i
  latlon[:,1]=np.arange(nsnodes)
  stations[vname]={"id":i,"nsnodes":nsnodes,"start":sIndex,"end":eIndex,"latlng":latlon,"featureid":np.zeros(nsnodes)+i}
  isnode=isnode+nsnodes