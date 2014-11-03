################################
# Job Generator pour MALBEC
# Il faut faire un script qui
# 1) Recupere les parametres a faire tourner
# 2) Creer une liste de dict pour remplir constante.f90
# 3) Boucle sur ces dicts:
#     - Creer un fichier constante temporaire
#     - Lancer l'execution du cript
#     - Creer un fichier job
#     - Lancer le job

################################
# Libray

import itertools
import subprocess
import sys
import os
from sys import platform as _platform
import distutils.util
import shutil
import datetime


################################
# 1) Parametre a creer :

if _platform == "linux" or _platform == "linux2":
    Root = '//gpfs/users/thorey/'
    Racine = '/home/thorey/Code_ELAS/'
    Root_Run = 'ELAS/Run_New_2/'
    Bactrack_Run = 'Bactrack_ELAS.txt'
    Compilateur = 'ifort'
elif _platform == "darwin":
    Root = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/'
    Racine = Root+'Code/'
    Root_Run = 'Code_ELAS/Test/Run/'
    Bactrack_Run = 'Bactrack_ELAS.txt'
    Compilateur = 'gfortran'
    
# Dict_Param = {'Sigma': ['2D-2'],
#               'Delta0': ['5D-3'],
#               'Grav': ['0D0'],
#               'El': ['1D0'],
#               'Nu': ['1D0','1D-2','1D-3'],
#               'Pe': ['1D0','1D-2','1D-3'],
#               'Psi': ['0.D0','3D-1'],
#               'N1' : ['1D5','1D0'],
#               'Dr' : ['1D-2'],
#               'Ep': ['1D-4'],
#               'Dt' : ['1D-6']}

Dict_Param = {'Sigma': ['2D-2'],
              'Delta0': ['5D-3'],
              'Grav': ['0D0'],
              'El': ['1D0'],
              'Nu': ['1D-2','1D-3'],
              'Pe': ['1D0','1D-1','1D-2'],
              'Psi': ['0.D0'],
              'N1' : ['1D5','1D0'],
              'Dr' : ['1D-2'],
              'Ep': ['1D-4'],
              'Dt' : ['1D-6']}


Init = 0 # 1 If you want to begin for the last backup
space = '\n --------------------- \n'

################################
# Date + fichier a ecrire + indications

now = datetime.datetime.now()
write = distutils.util.strtobool(input("Do you want to write in backtrace ? yes or no ?: "))
if write:
    with open(Racine + Bactrack_Run, 'a') as f:
        f.write('\n\n'+'-----------------------------'+'\n\n')
        f.write(str(now))
        f.write('\n\n'+'-----------------------------'+'\n\n')
        f.write('What are you running ?\n')
        f.write(str(raw_input('What are you running ?\n\n')) + '\n\n')
        f.write('Why do you want to test ?\n')
        f.write(str(raw_input('What do you want to test ? \n\n')) +'\n\n')
        
################################
# 2) Dictionnaire de dictionnaire

product = [x for x in apply(itertools.product, Dict_Param.values())]
Dict_Run = [dict(zip(Dict_Param.keys(), p)) for p in product]

################################
# 3) Boucle sur les runs


for run in Dict_Run:

    name = str('E' + run['El']
               + '_G' + run['Grav']
               + '_N' + run['Nu']
               + '_P' + run['Pe']
               + '_D' + run['Delta0']
               + '_C' + run['Psi']
               + '_R' + run['N1']
               + '_S' + run['Sigma']
               + '_Dr' + run['Dr']
               + '_Ep' + run['Ep']
               + '_Dt' + run['Dt'])

    print space + name + space

    if Init == 0:
        print ' Start from no backup \n'
        Backup = [ '0' , 'Backup_000000.dat']
        if os.path.isdir(Root+Root_Run+name):
            Bool =distutils.util.strtobool(input(" Directory already exist, Do you want to remove it ? yes or no ?: "))
            if Bool:
                shutil.rmtree(Root+Root_Run+name)
                os.mkdir(Root+Root_Run+name)
            else:
                sys.exit()
        else:
            os.mkdir(Root+Root_Run+name)
    elif Init == 1:
        if os.path.isdir(Root+Root_Run+name):
            if len([elt for elt in os.listdir(Root+Root_Run+name) if elt.split('_')[0] == 'Backup'])<2:
                print 'Pas de fichier backup dans le dossier '+name+' considerer. \n Consider to use Init=0'
                raise SystemExit
            else :
                back = [elt for elt in os.listdir(Root+Root_Run+name) if elt.split('_')[0] == 'Backup'][-1]
                Backup = [ '1' , back]
                print 'On repart du fichier ' + back 

        else:
            print 'Pas de fichie '+str(name)+' backup dans le dossier considerer, consider to use Init=0'
            raise SystemExit

    else :

        print 'Init can only take value between 1 ( Use of backup) and 0 ( from the beginning)'
        raise SystemExit
            
    # Creation de constante_tmp.f90 temporaire
    with open( str(Racine)+'Module_Init.f90' , 'r') as script:
            with open(str(Racine)+'Module_Init'+'_tmp.f90', 'wr+') as script_tmp:
                for l in script:
                    if l == '    CHARACTER(LEN=Size) :: Root\n':
                        to_write = l.replace('Size',str(len(Root)))
                    elif l == '    Sigma = Null\n':
                        to_write = l.replace('Null',run['Sigma'])
                    elif l == '    Delta0 = Null\n':
                        to_write = l.replace('Null',run['Delta0'])
                    elif l == '    Grav = Null\n':
                        to_write = l.replace('Null',run['Grav'])        
                    elif l == '    El = Null\n':
                        to_write = l.replace('Null',run['El'])
                    elif l == '    Nu = Null\n':
                        to_write = l.replace('Null',run['Nu'])
                    elif l == '    Pe = Null\n':
                        to_write = l.replace('Null',run['Pe'])
                    elif l == '    Psi = Null\n':
                        to_write = l.replace('Null',run['Psi'])
                    elif l == '    N1 = Null\n':
                        to_write = l.replace('Null',run['N1'])
                    elif l == '    Dr = Null\n':
                        to_write = l.replace('Null',run['Dr'])
                    elif l == '    eps_1 = Null\n':
                        to_write = l.replace('Null',run['Ep'])
                    elif l == '    Dt = Null\n':
                        to_write = l.replace('Null',run['Dt'])
                    elif l == '    Init = Null\n':
                        to_write = l.replace('Null',Backup[0])

                    elif l == '    Root = Null\n':
                        to_write = l.replace('Null', "'" + Root + "'")
                        
                    elif l == '    Input_Data_Name = Null\n':
                        to_write = l.replace('Null', "'" + Backup[1] + "'")
                        
                    elif l == '    Input_Racine = Root//Null\n':
                        to_write = l.replace('Null', "'" + Root_Run + name + "/'")
                        
                    elif l == '    Output_Racine = Root//Null\n':
                        to_write = l.replace('Null', "'" + Root_Run + name + "/'")
                        
                    else:
                        to_write = l
                    script_tmp.write(to_write)
                    
    # Make executble with name: name
    with open(str(Racine)+'run.sh' , 'r') as script:
            with open(str(Racine)+'run_tmp.sh', 'wr+') as script_tmp:
                for l in script:
                    if len(l.split('.o')) > 1:
                        to_write = l.replace('run',name)
                    else:
                        to_write = l
                    to_write = to_write.replace('Compilateur',Compilateur)
                    script_tmp.write(to_write)
                 
    subprocess.call(str(Racine)+'run_tmp.sh', shell=True)


    # Make the fihcier to run on malbec
    if _platform == "linux" or _platform == "linux2":
        with open( str(Racine) + 'run.job' , 'r') as script:
            with open(str(Racine)+name+'.job', 'wr+') as script_tmp:
                for l in script:
                    if l == '#SBATCH -J Run\n':
                        to_write = l.replace('Run',name)
                    elif l == './Run\n':
                        to_write = l.replace('Run',name)
                    elif l == 'cd dir\n':
                        to_write = l.replace('dir',Racine)
                    else:
                        to_write = l
                    script_tmp.write(to_write)
                        
        script = str('sbatch '+ Racine+ name+'.job')
        subprocess.call(script ,shell=True)

    if write:
        with open(Racine + Bactrack_Run, 'a') as f:
            f.write(name +'\n')
        
    print '\n  SUCESSFULL'
