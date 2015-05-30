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
    Root= '/gpfs/users/thorey/'
    Root_Code = '/home/thorey/EPGC/'
    Name_Folder_Run = '' # Remplir si on veut faire un test dans un dossier specific
    Bactrace_Run = 'Bactrack.txt'
    Compilateur = 'ifort'
elif _platform == "darwin":
    Root = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/TEST/'
    Root_Code = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/EPGC/'
    Name_Folder_Run = '' # Remplir si on veut faire un test dans un dossier specific
    Bactrace_Run = 'Bactrack.txt'
    Compilateur = 'gfortran'

# !Definition
# ! Model: {0: Integration Epaisseur, 1: Skin thermal layer}
# ! Schema :{0: Newton_Rhaspod, 1: Finite difference}
# ! Rheology: {0: Bercovici, 1: Roscoe, 2: Arrhenius}

Namejob = 'Arrhe_db2'
Model = 1
T_Schema = 0;H_Schema = 0
Rheology = 0

Dict_Param = {'Sigma': ['2D-2'],
              'Delta0': ['5D-3'],
              'Grav': ['0D0'],
              'El': ['1D0'],
              'Nu': ['1D0','1D-3'],
              'Pe': ['1D0','1D-3'],
              'Psi': ['1D0'],
              'N1' : ['1D0','1D5'],
              'gam':['0D0'],
              'Inter_Q':['1D20'],
              'Dr' : ['1D-2'],
              'Ep': ['1D-4'],
              'Dt' : ['1D-7']}
M_grid = 4000
Init = 0 # 1 If you want to begin for the last backup
space = '\n --------------------- \n'

def copy_folder(src,dest):
    os.mkdir(dest)
    print dest
    for filee in [f for f in os.listdir(src) if f[0]!='.']:
        print filee
        if not os.path.isdir(os.path.join(src,filee)):
            shutil.copy(filee,dest+'/')

def Destination_Folder(el,grav,Root):
    if el == 0.0 and grav == 0.0:
        print 'ERROR sur les parametres, el \
        and grav ne peuvent pas etre nulle en meme temps'
        sys.exit()
    elif el == '1D0' and grav == '0D0':
        return os.path.join(Root,'ELAS')
    elif el == '1D0' and grav == '1D0':
        return os.path.join(Root,'ELASGRAV')
    elif el == '0D0' and grav == '1D0':
        return os.path.join(Root,'GRAV')
    else:
        print 'el = %s, grav = %s'%(el,grav)
        print 'You messed the value of grav and el, do you?'
        sys.exit()

def Journal_name(el,grav):
    if el == '1D0' and grav == '0D0':
        return 'Journal_ELAS.txt'
    elif el == '1D0' and grav == '1D0':
        return 'Journal_ELAS_GRAV.txt'
    elif el == '0D0' and grav == '1D0':
        return 'Journal_GRAV.txt'
    else:
        print 'You messed the value of grav and el, do you?'
        sys.exit()
    
################################
# Creation d'un dossier pour acceuillier les simus
Dict_Model = {0: 'Int_Epaisseur', 1: 'Skin'}
Dict_Schema = {0: 'Newton', 1: 'GFD'}
Dict_Rheology = {0: 'Bercovici', 1: 'Roscoe', 2: 'Arrhenius'}
Root_Run = os.path.join(Destination_Folder(Dict_Param['El'][0],Dict_Param['Grav'][0],Root),'M'+Dict_Model[Model]+'_'\
                        'TSc_'+Dict_Schema[T_Schema]+'_'\
                        'HSc_'+Dict_Schema[H_Schema]+'_'\
                        'R'+Dict_Rheology[Rheology])
Journal = Journal_name(Dict_Param['El'][0],Dict_Param['Grav'][0])

if not os.path.isdir(Destination_Folder(Dict_Param['El'][0],Dict_Param['Grav'][0],Root)):
    os.mkdir(Destination_Folder(Dict_Param['El'][0],Dict_Param['Grav'][0],Root))
if not os.path.isdir(Root_Run):
    os.mkdir(Root_Run)

if Name_Folder_Run == '':
    print datetime.date.today()
    print 'On cree un nouveaux dossier pour acceuillir les runs, le code et le backtrace:'
    List_Folder_Compatible = [ f for f in os.listdir(Root_Run) if len(f.split('_')) == 3]
    List_Folder_Existant = [ f for f in List_Folder_Compatible  if f.split('_')[1] == str(datetime.date.today()) ]
    if len(List_Folder_Existant) == 0:
        Name_Folder_Run = 'Run_'+str(datetime.date.today())+'_0'
    else:
        Nombre_Folder = max([int(version.split('_')[2]) for version in List_Folder_Existant])
        Name_Folder_Run = 'Run_'+str(datetime.date.today())+'_'+str(Nombre_Folder+1)
    print 'Name_Folder_Run',Name_Folder_Run
    print 'Root_Code',Root_Code
    print os.path.join(Root_Run,Name_Folder_Run,'Run_Code')
    os.mkdir(os.path.join(Root_Run,Name_Folder_Run))
    copy_folder(Root_Code,os.path.join(Root_Run,Name_Folder_Run,'Run_Code'))
else :
    if not os.path.isdir(os.path.join(Root_Run,Name_Folder_Run)):
        print os.path.join(Root_Run,Name_Folder_Run)
        os.mkdir(os.path.join(Root_Run,Name_Folder_Run))
        copy_folder(os.path.join(Root_Code,Root_Run,Name_Folder_Run,'Run_Code'))
        # shutil.copytree(Root_Code,Root_Run+Name_Folder_Run+'/Run_Code')

Root_Run = os.path.join(Root_Run,Name_Folder_Run)
Root_Code = os.path.join(Root_Run,'Run_Code')

if not os.path.isdir(Root_Run) or not os.path.isdir(Root_Code):
    print 'Root_Run : ' + Root_Run
    print 'Root_Code : '+ Root_Code
    print 'Erreur dans un des deux paths ci dessus'
    print 'Vraissemblablement, il n existe pas !!!'
    sys.exit()
print 'Root_Run : ' + Root_Run
print 'Root_Code : ' + Root_Code

################################
# Journal - Record all the version runed
if Init == 0 :
    now = datetime.datetime.now()
    write = distutils.util.strtobool(input("Do you want to write in the journal ?: "))
    if write:
        with open(os.path.join(Destination_Folder(Dict_Param['El'][0],Dict_Param['Grav'][0],Root),Journal), 'a') as f:
            f.write('\n'+'####################'+'\n')
            f.write('-----------------------------'+'\n')
            f.write(str(now)+' ------- '+Name_Folder_Run)
            f.write('\n'+'-----------------------------'+'\n')
            f.write('Short Description of the runs you are aiming to run ?\n')
            f.write(str(raw_input('Short descriptin of '+Name_Folder_Run+' ?\n'+'\n')))
            f.write('\n'+'-----------------------------'+'\n')
            f.write('Model = '+Dict_Model[Model]+'\n')
            f.write('T_Schema = '+Dict_Schema[T_Schema]+'\n')
            f.write('H_Schema = '+Dict_Schema[H_Schema]+'\n')
            f.write('Rheology = '+Dict_Rheology[Rheology]+'\n')
            f.write('\n'+'-----------------------------'+'\n')
            f.write('Liste des parametres, i.e. Nombre de Runs'+'\n')
            for key,item in Dict_Param.iteritems():
                f.write(key+': '+str(item)+'\n')
            f.write('\n'+'####################'+'\n')
            
################################
# Date + fichier a ecrire + indications
    
now = datetime.datetime.now()
write = distutils.util.strtobool(input("Do you want to write in backtrace ? yes or no ?: "))
if write:
    with open(os.path.join(Root_Code,Bactrace_Run), 'a') as f:
        f.write('\n\n'+'-----------------------------'+'\n\n')
        f.write(str(now))
        f.write('\n\n'+'-----------------------------'+'\n\n')
        f.write('What and Why are you running ?\n')
        f.write(str(raw_input('What are you running ?\n\n')) + '\n\n')
        f.write('Runs lancer sur '+_platform +'\n\n')
        
################################
# 2) Dictionnaire de dictionnaire

product = [x for x in apply(itertools.product, Dict_Param.values())]
Dict_Run = [dict(zip(Dict_Param.keys(), p)) for p in product]
print 'On s"apprete a lancer '+str(len(Dict_Run))+' jobs'

################################
# 3) Boucle sur les runs
compteur_run = len(Dict_Run)
for run in Dict_Run:
    print 'Il reste rencore %d runs a lancer'%(compteur_run)
    compteur_run -= 1
    name = str('E' + run['El']
               + '_G' + run['Grav']
               + '_N' + run['Nu']
               + '_P' + run['Pe']
               + '_D' + run['Delta0']
               + '_C' + run['Psi']
               + '_R' + run['N1']
               + '_G' + run['gam']
               + '_I' + run['Inter_Q']
               + '_S' + run['Sigma']
               + '_Dr' + run['Dr']
               + '_Ep' + run['Ep']
               + '_Dt' + run['Dt'])

    print space + name + space

    Backup = [ '0' , 'Backup_000000.dat']
    if os.path.isdir(os.path.join(Root_Run,name)):
        print 'Le repertoire existe deja: Voici la liste des fichiers:'
        print [elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup']
        if len([elt for elt in os.listdir(os.path.join(Root_Run,name)) if elt.split('_')[0] == 'Backup']) < 10:
            Bool = True
        else:
            Bool = False
            # Bool =distutils.util.strtobool(input(" Directory already exist, Do you want to remove it ? yes or no ?: ")) 
        if Bool:
            shutil.rmtree(os.path.join(Root_Run,name))
            os.mkdir(os.path.join(Root_Run,name))
        else:
            if len([elt for elt in os.listdir(os.path.join(Root_Run,name)) if elt.split('_')[0] == 'Backup'])<2:
                print 'Pas de fichier backup dans le dossier '+name+' considerer. \n Consider to say do to the question before'
                raise SystemExit
            else :
                back = [elt for elt in os.listdir(os.path.join(Root_Run,name)) if elt.split('_')[0] == 'Backup'][-1]
                Backup = [ '1' , back]
                print 'On repart du fichier ' + back 
    else:
        os.mkdir(os.path.join(Root_Run,name))

            
    # Creation de constante_tmp.f90 temporaire
    with open( os.path.join(Root_Code,'Module_Init.f90') , 'r') as script:
            with open(os.path.join(Root_Code,'Module_Init'+'_tmp_1.f90'), 'wr+') as script_tmp:
                for l in script:
                    if l == '    CHARACTER(LEN=Size) :: Root\n':
                        to_write = l.replace('Size',str(len(Root_Run)))
                    elif l == '    CHARACTER(LEN=Size_Name) :: Name_File\n':
                        to_write = l.replace('Size_Name',str(len(name)))
                    elif l == '    CHARACTER(LEN=Size) :: Root_Code\n':
                        to_write = l.replace('Size',str(len(Racine)))
                    elif l == '    Model = Null\n':
                        to_write = l.replace('Null',str(Model))
                    elif l == '    T_Schema = Null\n':
                        to_write = l.replace('Null',str(T_Schema))
                    elif l == '    H_Schema = Null\n':
                        to_write = l.replace('Null',str(H_Schema))
                    elif l == '    Rheology = Null\n':
                        to_write = l.replace('Null',str(Rheology))
                    elif l == '    M = Null\n':
                        to_write = l.replace('Null',str(M_grid))
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
                    elif l == '    gam = Null\n':
                        to_write = l.replace('Null',run['gam'])
                    elif l == '    Inter_Q = Null\n':
                        to_write = l.replace('Null',run['Inter_Q'])
                    elif l == '    Dr = Null\n':
                        to_write = l.replace('Null',run['Dr'])
                    elif l == '    eps_1 = Null\n':
                        to_write = l.replace('Null',run['Ep'])
                    elif l == '    Dt = Null\n':
                        to_write = l.replace('Null',run['Dt'])
                    elif l == '    NF = Null\n':
                        to_write = l.replace('Null',"'"+name+"'")
                    elif l == '    Init = Null\n':
                        to_write = l.replace('Null',Backup[0])
                    elif l == '    Root = Null\n':
                        to_write = l.replace('Null', "'" + Root_Run + "/'")
                    elif l == '    Root_Code = Null\n':
                        to_write = l.replace('Null', "'" + Root_Code + "/'")
                        
                    elif l == '    Input_Data_Name = Null\n':
                        to_write = l.replace('Null', "'" + Backup[1] + "'")
                        
                    elif l == '    Input_Racine = Root//Null\n':
                        to_write = l.replace('Null', "'/" +  name + "/'")
                        
                    elif l == '    Output_Racine = Root//Null\n':
                        to_write = l.replace('Null', "'/" + name + "/'")
                    else:
                        to_write = l
                    script_tmp.write(to_write)

    # ON decoupe les lignes trop grandes
    with open( os.path.join(Root_Code,'Module_Init_tmp_1.f90'), 'r') as script:
        with open(os.path.join(Root_Code,'Module_Init'+'_tmp.f90'), 'wr+') as script_tmp:
            for line in script:
                if line.split('=')[0] != '    Root ' and line.split('=')[0] != '    Root_Code ':
                    to_write = line
                    script_tmp.write(to_write)
                else:
                    to_write1 = line[:70]+"'&\n"
                    to_write2 = "&//'"+line[70:]
                    script_tmp.write(to_write1)
                    script_tmp.write(to_write2)
            
    # Make executble with name: name
    with open(os.path.join(Root_Code,'run.sh') , 'r') as script:
            with open(os.path.join(Root_Code,'run_tmp.sh'), 'wr+') as script_tmp:
                for l in script:
                    if len(l.split('.o')) > 1:
                        to_write = l.replace('run',os.path.join(Root_Code,name))
                    else:
                        to_write = l
                    to_write = to_write.replace('Compilateur',Compilateur)
                    script_tmp.write(to_write)

    subprocess.call(os.path.join(Root_Code,'run_tmp.sh'), shell=True)

# Make the fihcier to run on malbec
if _platform == "linux" or _platform == "linux2":
    with open( os.path.join(Root_Code, 'run.job') , 'r') as script:
        with open(os.path.join(Root_Code,'run_tmp'+'.job'), 'wr+') as script_tmp:
            for l in script:
                if l == '#SBATCH -J Name\n':
                    to_write = l.replace('Name',Namejob)
                elif l == '#SBATCH --nodes=Nnode\n':
                    to_write = l.replace('Nnode',str((len(Dict_Run)-1)//16+1))
                elif l == '#SBATCH --ntasks=Ntask\n':
                    to_write = l.replace('Ntask',str((len(Dict_Run))))
                elif l == 'cd dir\n':
                    to_write = l.replace('dir',Root_Code)
                else:
                    to_write = l
                script_tmp.write(to_write)
    #On ecrit le fichier Multi_Job
    map(os.remove,[os.path.join(Root_Code,f) for f in os.listdir(str(Root_Code))
                   if f.split('_')[0] == 'Multi'])
    Job_File = os.path.join(Root_Code,'Multi_Run')
    for i,run in enumerate(Dict_Run):
        name = str('E' + run['El']
                   + '_G' + run['Grav']
                   + '_N' + run['Nu']
                   + '_P' + run['Pe']
                   + '_D' + run['Delta0']
                   + '_C' + run['Psi']
                   + '_R' + run['N1']
                   + '_G' + run['gam']
                   + '_I' + run['Inter_Q']
                   + '_S' + run['Sigma']
                   + '_Dr' + run['Dr']
                   + '_Ep' + run['Ep']
                   + '_Dt' + run['Dt'])
        with open(Job_File, 'a') as fc:
            fc.write(str(i)+' '+name+'\n')
            with open(os.path.join(Root_Code ,Bactrace_Run), 'a') as f:
                f.write(name +'\n')

    script = str('sbatch '+ os.path.join(Root_Code,'run_tmp'+'.job'))
    subprocess.call(script ,shell=True)
