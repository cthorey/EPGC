##############
# Library

from sys import platform as _platform
import os , sys , subprocess , itertools
import cPickle as pickle
import pandas as pd
import shutil
import numpy as np
from itertools import groupby
##############
# Fonction

def set_path_input(Root,Run):
    run = os.path.join(Root,Run)
    if not os.path.exists(run):
        print 'Erreur, pas de run en question',run
        sys.exit()
    return run
        
def set_path_output(Root,Run):
    workspace = (os.path.join(Root,Run)).replace('Run','Workspace')
    if not os.path.exists(workspace):
        os.mkdir(workspace)
    return workspace
    
def Who_is_Root():
    if _platform == "linux" or _platform == "linux2":
        Root = '/gpfs/users/thorey/'
    elif _platform == "darwin":
        Root = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/'

    return Root

def to_number(s):
    try:
        s1 = np.float(s)
        return s1
    except ValueError:
        return np.float(0)

def m_Test_Size_pickle(input_file):
        with open(input_file, 'rb') as f:
            try:
                data_pickle = pd.DataFrame(pickle.load(f)['Max'])
            except:
                data_pickle = pd.DataFrame(pd.read_pickle(input_file)['Max'])
        return len(data_pickle.tm)
        
def Tcheque_Workspace(key,Workspace,c_input,c_output):
    filee = os.path.join(c_output,Workspace)
    boole= None
    if not os.path.isfile(filee):
        # print 'N existe pas encore',filee
        boole = False
    else:
        N = m_Test_Size_pickle(filee)
        N_new = len([f for f in os.listdir(os.path.join(c_input,key)) if f.split()[0][0] == 'R'])
        # print N,N_new
        if N == N_new:
            boole = True
        else:
            boole = False
    return boole
    
def to_string(s):
    return str(s).replace('D','E')

def Load_Nsd(Racine):
    test = os.path.isfile(os.path.join(Racine,'NbSsDim.txt'))
    if not test:
        # print 'Pas de fichier: '+ key
        return test,None,None
    Nsd = pd.read_csv(os.path.join(Racine,'NbSsDim.txt'),
                      sep=' ',
                      skipinitialspace=True,
                      lineterminator='\n',
                      header=0)
    Nsd = Nsd.applymap(lambda x: to_string(x))
    Nsd = Nsd.applymap(lambda x: to_number(x))
    if float(Nsd.Pe) == 0:
        Nsd['Pe'] = 0.0
    else:
        Nsd['Pe']= 1./float(Nsd.Pe)
    Workspace = '_'.join([f[:2]+str(np.float(Nsd[f][0])) for f in Nsd.columns.tolist()])

    return test,Workspace,Nsd

def load_Data_tmp(Racine):
    Data = []
    files = [f for f in os.listdir(Racine) if f.split()[0][0] == 'R']
    for i, file in enumerate(files):
        Data_tmp = pd.read_csv(os.path.join(Racine,file),
                               sep=' ',
                               skipinitialspace=True,
                               lineterminator='\n',
                               header=0)
        if len(Data_tmp) == 0:
            os.remove(os.path.join(Racine,file))
            os.remove(os.path.join(Racine,'Backup_'+file[3:]))
            continue
        Data_tmp = Data_tmp.applymap(lambda x: to_string(x))
        Data_tmp = Data_tmp.applymap(lambda x: to_number(x))
        Data.append(Data_tmp)

    return Data

def convert_to_worskpace(c_input,c_output,tracker):
    subdirs = [f for f in os.listdir(c_input)
               if (os.path.isdir(os.path.join(c_input, f)) and f[0] == 'E')]
    merge_d = {f[:-7]:f for f in subdirs}.keys()
    merge = [[g for g in subdirs if g[:-7] == f ] for f in merge_d]

    track = open(tracker,'w')
    track.write('Dossier input %s\n'%(c_input))
    track.close()
    for liste in merge:
        print liste
        try:
            track = open(tracker,'a')
            track.write(liste[0]+'\n')
            track.close()
            dict_tmp = dict.fromkeys(liste)
            if len(dict_tmp) == 1: # On fait un premier test pour savoir si le wroksapce existe deja
                key = dict_tmp.keys()[0]
                Racine = os.path.join(c_input,key)
                test,Workspace,Nsd = Load_Nsd(Racine)
                Already_Same_Workspace = Tcheque_Workspace(key,Workspace,c_input,c_output)
                if Already_Same_Workspace:
                    # print 'Pas de changement dans: '+key
                    continue
            for key,value in dict_tmp.iteritems():
                Racine = os.path.join(c_input,key)
                test,Workspace,Nsd = Load_Nsd(Racine)
                Data = load_Data_tmp(Racine)
                dict_tmp[key] = {'NsD': Nsd,
                                 'Data': Data}
        
            data_tmp = []; Nsd_l = [];
            for df in dict_tmp.itervalues():
                Nsd_l.append(df['NsD'])
                data_tmp.append(df['Data'])
            Data_final=[df for df in itertools.chain(*data_tmp)]
            t_index=np.array([[i,df.tm.max()] for i,df in enumerate(Data_final)])
            df=pd.DataFrame({'indice':t_index[:,0],'t':t_index[:,1]}).sort('t')
            indice=[np.int(elt) for elt in df.sort('t').indice]
            Data_final=[Data_final[i] for i in indice]
            Column_data = ['tm','dist', 'H','Te','BL', 'Xi','Ts', 'P','Srr','Stt','R','hmubar','hthetabar','Mu_e']
            Data = [Df[Column_data] for Df in Data_final]
            for i,df in enumerate(Data_final):
                if i == 0:
                    Max = { key : [] for key in df.columns}
                for key,liste in Max.iteritems():
                    liste.append(df[key][0])
            
            D_pickle = {'NsD': Nsd_l[0],
                        'Data': Data,
                        'Max': Max}
            with open(os.path.join(c_output,Workspace), 'wb') as f:
                pickle.dump(D_pickle, f, pickle.HIGHEST_PROTOCOL)
        except:
            track = open(tracker,'a')
            track.write(liste[0]+': FAILED \n')
            track.close()        

def copy_folder(src,dest):
    if not os.path.isdir(dest):
        os.mkdir(dest)
    for filee in [f for f in os.listdir(src) if f[0]!='.']:
        if not os.path.isdir(os.path.join(src,filee)):
            shutil.copy(os.path.join(src,filee),dest)                

#############
# Program

root_path = Who_is_Root()
runs = []
# runs.append(os.path.join('ELAS','MSkin_TSc_Newton_HSc_Newton_RArrhenius','Run_2015-05-30_1'))
runs.append(os.path.join('ELAS','MSkin_TSc_Newton_HSc_Newton_RArrhenius','Run_2015-05-30_4'))
# runs.append(os.path.join('ELAS','MSkin_TSc_Newton_HSc_Newton_RBercovici','Run_2015-05-30_1'))
# runs.append(os.path.join('GRAV','MSkin_TSc_GFD_HSc_Newton_RArrhenius','Run_2015-05-30_0'))
# runs.append(os.path.join('GRAV','MSkin_TSc_GFD_HSc_Newton_RBercovici','Run_2015-05-30_0'))

# runs.append('GRAV/Run_2015-04-21_0/')# Pas obulier / a la fin
for run in runs:
    tracker = '/home/thorey/EPGC/Tracker_'+'_'.join(run.split('/'))+'.txt'
    run_path = set_path_input(root_path,run)
    workspace_path = set_path_output(root_path,run)
    copy_folder(os.path.join(run_path,'Run_Code'),os.path.join(workspace_path,'Run_Code'))
    convert_to_worskpace(run_path,workspace_path,tracker)
    f = open(tracker,'a')
    f.write('Workspace made with sucess !! \n')
    f.close()
    
    
