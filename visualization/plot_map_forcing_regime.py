import sys
import os
sys.path.append(os.environ['GOTMWORK_ROOT']+'/tools', )
from gotmanalysis import *
np.seterr(all='raise')

# timetag = '20080701-20080731'
timetag = '20090101-20090131'
# casename = 'COREII_Global'
casename = 'JRA55-do_Global'
forcing_reg_type = 'BG12'
# forcing_reg_type = 'LF17'

update_data = False

s1data_root = os.environ['GOTMRUN_ROOT']+'/'+casename+'/VR1m_DT600s_'+timetag
s2data_root = os.environ['GOTMFIG_ROOT']+'/data/'+casename+'/VR1m_DT600s_'+timetag
fig_root = os.environ['GOTMFIG_ROOT']+'/'+casename+'/VR1m_DT600s_'+timetag
os.makedirs(s2data_root, exist_ok=True)
os.makedirs(fig_root, exist_ok=True)

tmname = 'KPP-CVMix'
basepath = s1data_root+'/'+tmname
s2data_name = s2data_root+'/data_forcing_regime_'+forcing_reg_type+'_'+tmname+'.npz'
figname = fig_root+'/fig_forcing_regime_'+forcing_reg_type+'.png'
loclist = sorted(os.listdir(basepath))
if update_data or not os.path.isfile(s2data_name):
    # save data
    pathlist = [basepath+'/'+x+'/gotm_out_s1.nc' for x in loclist]
    godmobj = GOTMOutputDataMap(pathlist)
    forcing_regime = np.zeros(godmobj.ncase)
    unstable = np.zeros(godmobj.ncase)
    for i in np.arange(godmobj.ncase):
        if np.mod(i, 100) == 0:
            print('{:6.2f} %'.format(i/godmobj.ncase*100.0))
        tmp = GOTMOutputData(godmobj._paths[i], init_time_location=False)
        if forcing_reg_type == 'BG12':
            forcing_regime[i] = tmp.diag_forcing_regime_BG12()
        elif forcing_reg_type == 'LF17':
            forcing_regime[i] = tmp.diag_forcing_regime_LF17()

    gmobj = GOTMMap(data=forcing_regime, lon=godmobj.lon, lat=godmobj.lat, name='forcing_regime')
    gmobj.save(s2data_name)
else:
    # read data
    gmobj = GOTMMap().load(s2data_name)


fig = plt.figure()
fig.set_size_inches(6, 2.2)
plot_forcing_regime(gmobj)
plt.tight_layout()
plt.savefig(figname, dpi = 300)

