# shared variables for COREII and JRA55do cases

# list of turbulence closure methods (for reference to the paths)
list_turbmethod = ['KPP-CVMix',
                   'KPP-ROMS',
                   'KPPLT-EFACTOR',
                   'KPPLT-ENTR',
                   'KPPLT-RWHGK',
                   'EPBL-RH18',
                   'EPBL-RL19',
                   'SMC',
                   'SMCLT',
                   'K-EPSILON-SG',
                   'SMC-C01A',
                   'OSMOSIS']

# list of name of the turbulence closure methods (for display)
list_tmname = ['KPP-CVMix',
               'KPP-ROMS',
               'KPPLT-VR12',
               'KPPLT-LF17',
               'KPPLT-R16',
               'ePBL',
               'ePBL-LT',
               'SMC-KC94',
               'SMCLT-H15',
               '$k$-$\epsilon$-SG95',
               'SMC-C01A',
               'OSMOSIS']

# list of simulation setup
list_setup = ['VR1m_DT600s_20090101-20090131',
              'VR1m_DT600s_20090201-20090228',
              'VR1m_DT600s_20090301-20090331',
              'VR1m_DT600s_20090401-20090430',
              'VR1m_DT600s_20090501-20090531',
              'VR1m_DT600s_20080601-20080630',
              'VR1m_DT600s_20080701-20080731',
              'VR1m_DT600s_20080801-20080831',
              'VR1m_DT600s_20080901-20080930',
              'VR1m_DT600s_20081001-20081031',
              'VR1m_DT600s_20081101-20081130',
              'VR1m_DT600s_20081201-20081231']

# list of month names
list_month = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']

# for plotting maps in two columns
irow_2col = [1, 2, 0, 1, 2, 3, 3, 4, 4, 5, 5]
icol_2col = [0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1]
labels_2col = ['(b)', '(c)', '(g)', '(h)', '(i)', '(d)', '(j)', '(e)', '(k)','(f)','(l)']
l_nlt = [True, True, False, False, False, True, False, True, False, True, False]
