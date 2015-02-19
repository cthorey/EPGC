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
    Root_ELAS = '/gpfs/users/thorey/ELAS/'
    Root_Run = '/gpfs/users/thorey/ELAS/' # Modifier la pour le dossier ou l'on va travailler
    Root_Code = '/home/thorey/Code_ELAS/'
    Name_Folder_Run = '' # Remplir si on veut faire un test dans un dossier specific
    Bactrace_Run = 'Bactrack.txt'
    Journal_ELAS = 'Journal_ELAS.txt'
    Compilateur = 'ifort'
elif _platform == "darwin":
    Root_ELAS = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/Code_ELAS/TEST/Run/'
    Root_Run = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/Code_ELAS/TEST/Run/'
    Root_Code = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/Code_ELAS/'
    Name_Folder_Run = '' # Remplir si on veut faire un test dans un dossier specific
    Bactrace_Run = 'Bactrack.txt'
    Journal_ELAS = 'Journal_ELAS.txt'
    Compilateur = 'gfortran'

Dict_Param = {'Sigma': ['5D-2'],
              'Delta0': ['5D-3'],
              'Grav': ['0D0'],
              'El': ['1D0'],
              'Nu': ['1D0'],
              'Pe': ['0D0'],
              'Psi': ['0.D0','1D0','5D0'],
              'N1' : ['1D5'],
              'Dr' : ['1D-2'],
              'Ep': ['1D-4'],
              'Dt' : ['1D-6']}

M_grid = 4000
Init = 0 # 1 If you want to begin for the last backup
space = '\n --------------------- \n'

def copy_folder(src,dest):
    os.mkdir(dest)
    print dest
    for filee in [f for f in os.listdir(src) if f[0]!='.']:
        shutil.copy(filee,dest+'/')
        
################################
# Creation d'un dossier pour acceuillier les simus
if Init == 0:
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
        os.mkdir(Root_Run+Name_Folder_Run)
        copy_folder(Root_Code,Root_Run+Name_Folder_Run+'/Run_Code')
        # shutil.copytree(Root_Code,Root_Run+Name_Folder_Run+'/Run_Code')
    else :
        if not os.path.isdir(Root_Run+Name_Folder_Run):
            os.mkdir(Root_Run+Name_Folder_Run)
            shutil.copytree(Root_Code,Root_Run+Name_Folder_Run+'/Run_Code')
        print Root_Run+Name_Folder_Run
elif Init == 1:
    print Root_Run+Name_Folder_Run
    if Name_Folder_Run == '':
        print 'Preciser le Name_Folder_Run dans lequelle il faut travaille'
        print 'Init =1 implique que des donner ont deja etait produite'
        sys.exit()
else:
    print 'Init ne peut etre que 0 ou 1'
    sys.exit()
Root_Run = Root_Run+Name_Folder_Run+'/'
Root_Code = Root_Run+'Run_Code/'

if not os.path.isdir(Root_Run) or not os.path.isdir(Root_Code):
    print 'Root_Run : ' + Root_Run
    print 'Root_Code : '+ Root_Code
    print 'Erreur dans un des deux paths ci dessus'
    print 'Vraissemblablement, il n existe pas !!!'
    sys.exit()
print 'Root_Run : ' + Root_Run
print 'Root_Code : ' + Root_Code

################################
# Journal de ELAS- Record all the version runed
if Init == 0 :
    now = datetime.datetime.now()
    write = distutils.util.strtobool(input("Do you want to write in the journal ?: "))
    if write:
        with open(Root_ELAS+Journal_ELAS, 'a') as f:
            f.write('\n'+'####################'+'\n')
            f.write('-----------------------------'+'\n')
            f.write(str(now)+' ------- '+Name_Folder_Run)
            f.write('\n'+'-----------------------------'+'\n')
            f.write('Short Description of the runs you are aiming to run ?\n')
            f.write(str(raw_input('Short descriptin of '+Name_Folder_Run+' ?\n'+'\n')))
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
    with open(Root_Code+Bactrace_Run, 'a') as f:
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
        if os.path.isdir(Root_Run+name):
            print 'Le repertoire existe deja: Voici la liste des fichiers:'
            print [elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup']
            if len([elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup']) == 0:
                Bool = True
            else:
                Bool =distutils.util.strtobool(input(" Directory already exist, Do you want to remove it ? yes or no ?: ")) 
            if Bool:
                shutil.rmtree(Root_Run+name)
                os.mkdir(Root_Run+name)
            else:
                if len([elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup'])<2:
                    print 'Pas de fichier backup dans le dossier '+name+' considerer. \n Consider to say do to the question before'
                    raise SystemExit
                else :
                    back = [elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup'][-1]
                    Backup = [ '1' , back]
                    print 'On repart du fichier ' + back 
        else:
            os.mkdir(Root_Run+name)
    elif Init == 1:
        if os.path.isdir(Root_Run+name):
            if len([elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup'])<2:
                print 'Pas de fichier backup dans le dossier '+name+' considerer. \n Consider to use Init=0'
                raise SystemExit
            else :
                back = [elt for elt in os.listdir(Root_Run+name) if elt.split('_')[0] == 'Backup'][-1]
                Backup = [ '1' , back]
                print 'On repart du fichier ' + back 
        else:
            print 'Pas de fichie '+str(name)+' backup dans le dossier considerer, consider to use Init=0'
            raise SystemExit

    else :
        print 'Init can only take value between 1 ( Use of backup) and 0 ( from the beginning)'
        raise SystemExit
            
    # Creation de constante_tmp.f90 temporaire
    with open( str(Root_Code)+'Module_Init.f90' , 'r') as script:
            with open(str(Root_Code)+'Module_Init'+'_tmp.f90', 'wr+') as script_tmp:
                for l in script:
                    if l == '    CHARACTER(LEN=Size) :: Root\n':
                        to_write = l.replace('Size',str(len(Root_Run)))
                    elif l == '    CHARACTER(LEN=Size_Name) :: Name_File\n':
                        to_write = l.replace('Size_Name',str(len(name)))
                    elif l == '    CHARACTER(LEN=Size) :: Root_Code\n':
                        to_write = l.replace('Size',str(len(Racine)))
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
                        to_write = l.replace('Null', "'" + Root_Run + "'")
                    elif l == '    Root_Code = Null\n':
                        to_write = l.replace('Null', "'" + Root_Code + "'")
                        
                    elif l == '    Input_Data_Name = Null\n':
                        to_write = l.replace('Null', "'" + Backup[1] + "'")
                        
                    elif l == '    Input_Racine = Root//Null\n':
                        to_write = l.replace('Null', "'" +  name + "/'")
                        
                    elif l == '    Output_Racine = Root//Null\n':
                        to_write = l.replace('Null', "'" + name + "/'")
                    else:
                        to_write = l
                    script_tmp.write(to_write)

    # Make executble with name: name
    with open(str(Root_Code)+'run.sh' , 'r') as script:
            with open(str(Root_Code)+'run_tmp.sh', 'wr+') as script_tmp:
                for l in script:
                    if len(l.split('.o')) > 1:
                        to_write = l.replace('run',Root_Code+name)
                    else:
                        to_write = l
                    to_write = to_write.replace('Compilateur',Compilateur)
                    script_tmp.write(to_write)

    subprocess.call(str(Root_Code)+'run_tmp.sh', shell=True)

# Make the fihcier to run on malbec
if _platform == "linux" or _platform == "linux2":
    with open( str(Root_Code) + 'run.job' , 'r') as script:
        with open(str(Root_Code)+'run_tmp'+'.job', 'wr+') as script_tmp:
            for l in script:
                if l == '#SBATCH -J Name\n':
                    to_write = l.replace('Name','EL')
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
    map(os.remove,[str(Root_Code)+f for f in os.listdir(str(Root_Code))
                   if f.split('_')[0] == 'Multi'])
    Job_File = str(Root_Code)+'Multi_Run'
    for i,run in enumerate(Dict_Run):
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
        with open(Job_File, 'a') as fc:
            fc.write(str(i)+' '+name+'\n')
            with open(Root_Code + Bactrace_Run, 'a') as f:
                f.write(name +'\n')

    script = str('sbatch '+ str(Root_Code)+'run_tmp'+'.job')
    subprocess.call(script ,shell=True)
