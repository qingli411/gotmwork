# remove runs with NaNs in the output data

import sys
import os
import shutil
sys.path.append(os.environ['GOTMWORK_ROOT']+'/tools', )
from gotmanalysis import *
np.seterr(all='raise')

timetag = '20080701-20080731'
# timetag = '20090101-20090131'
casename = 'COREII_Global'
# casename = 'JRA55-do_Global'
s1data_root = os.environ['GOTMRUN_ROOT']+'/'+casename+'/VR1m_DT600s_'+timetag

tmname = 'KPP-CVMix'
basepath = s1data_root+'/'+tmname
loclist = sorted(os.listdir(basepath))
# save data
pathlist = [basepath+'/'+x+'/gotm_out_s1.nc' for x in loclist]
godmobj = GOTMOutputDataMap(pathlist)
filter_list = []
for i in np.arange(godmobj.ncase):
    if np.mod(i, 100) == 0:
        print('{:6.2f} %'.format(i/godmobj.ncase*100.0))
    tmp = GOTMOutputData(godmobj._paths[i], init_time_location=False)
    try:
        tmp.open()
        rho = tmp.dataset.variables['rho'][:]
        tmp.close()
    except FloatingPointError:
        filter_list.append(i)
print('Done')

print('Removing the following directories...')
tmlist = sorted(os.listdir(s1data_root))
ntm = len(tmlist)
for i in filter_list:
    locname = loclist[i]
    for tm in tmlist:
        filter_dir = s1data_root+'/'+tm+'/'+locname
        print(filter_dir)
        shutil.rmtree(filter_dir)

