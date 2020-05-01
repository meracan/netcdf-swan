from .dataTest import elem,time,lat,lon,bed,slat,slon,freq,dir,sgroup,spcgroup


def create_mat(filePath,dic):
  # dic:{"datetime":np.array(ntime,datetime64),"hs":np.array((ntime,nnode))}
  None
  
def create_spc(filePath,dic):
  # dic:{"datetime":datetime64,"datetime":np.array((ntime,nsnode,nfreq,ndir))}
  None

def create_data():
  # for year in year
    # for month in months:
      # for station in stations:
        # nodeperstation=stations[station]
        # create_spc      
      # for var in sgroup:
        # create_mat
        