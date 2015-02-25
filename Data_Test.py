#LIBRAIRIE A IMPORTER

import os
import pandas as pd
import numpy as np
import sys
import cPickle as pickle

Base = '/Users/thorey/Documents/These/Projet/Refroidissement/'
Base += 'Skin_Model/Code/Code_ELAS/'
# Base = '/Users/thorey/Documents/These/Projet/Refroidissement/modelisation/Current_Version/'
# Base += 'Temperature_Boundary_Layer/Code/Old_Version/Version_SCAPAD/'
c_input = Base + 'TEST/Run/'
c_output  = Base + 'TEST/Workspace/'
#c_input = Base + 'SCAPAD/GRAV/Run/'
#c_output = Base + 'SCAPAD/GRAV/Workspace/'
select = ['G1D0']


def to_number(s):
    try:
        s1 = np.float(s)
        return s1
    except ValueError:
        return np.float(0)

def to_string(s):
    return str(s).replace('D','E')

subdirs=['1']

for elt in subdirs:
    print elt
    Racine = c_input
    files = [f for f in os.listdir(Racine) if f.split()[0][0] == 'R']
    N = len(files)-1
    Data = []
    
    D_pickle = {}
    Column_Nsd = ['el', 'grav', 'delta0', 'sigma', 'nu', 'Pe','Psi',
                  'N1','M','Dt','Dr','eps']
    test = os.path.isfile(Racine+'NbSsDim.txt')
    if not test:
        continue
    Nsd = pd.read_csv(Racine+'NbSsDim.txt',
                      sep=' ',
                      skipinitialspace=True,
                      lineterminator='\n',
                      header=0)
    Nsd = Nsd.applymap(lambda x: to_string(x))
    Nsd = Nsd.applymap(lambda x: to_number(x))

    for i, file in enumerate(files):
        print file
        Data_tmp = pd.read_csv(Racine+file,
                               sep=' ',
                               skipinitialspace=True,
                               lineterminator='\n',
                               header=0)
        Data_tmp = Data_tmp.applymap(lambda x: to_string(x))
        Data_tmp = Data_tmp.applymap(lambda x: to_number(x))
        if i == 0:
            column_data = ['tm','dist', 'H','Te','BL', 'Xi','Ts', 'P','Srr','Stt','R','hmubar','hthetabar','Mu_e','ubar']
            Max = { key : [] for key in Data_tmp.columns}
        Data.append(Data_tmp[column_data])
        for key,liste in Max.iteritems():
            liste.append(Data_tmp[key][0])
            
    D_pickle = {'NsD': Nsd,
                'Data': Data,
                'Max': Max}
        
    Workspace = str('E' + str(np.float(Nsd.el[0]))
                    + '_G' + str(np.float(Nsd.grav[0]))
                    + '_N' + str(np.float(Nsd.nu[0]))
                    + '_P' + str(np.float(Nsd.Pe[0]))
                    + '_D' + str(np.float(Nsd.delta0[0]))
                    + '_C' + str(np.float(Nsd.Psi[0]))
                    + '_R' + str(np.float(Nsd.N1[0]))
                    + '_S' + str(np.float(Nsd.sigma[0]))
                    + '_Dr' + str(np.float(Nsd.Dr[0]))
                    + '_Ep' + str(np.float(Nsd.eps[0])))

    with open(c_output+Workspace, 'wb') as f:
        pickle.dump(D_pickle, f, pickle.HIGHEST_PROTOCOL)

