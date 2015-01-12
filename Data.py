
# Recupere les donne sur malbec,
# Load les donnes

import subprocess
import sys
import os

reset = 0
Folder = '2015-01-09_0'
Folder_malbec = 'Run_'+Folder+'/'
Folder_laptop = 'Run_'+Folder+'/'
Workspace_laptop = 'Workspace_'+Folder+'/'
Routine_python = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Routine/ELAS/E_Load'
Root_malbec = '/gpfs/users/thorey/ELAS/'+Folder_malbec
Root_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/ELAS/'+Folder_laptop
Workspace_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/ELAS/'+Workspace_laptop

# On fait un peu de menage
if reset == 1: # On enleve le dossier et tous ls workspace pr recommecner de zero
    try :
        subprocess.call('rm -R '+Root_laptop +'*',shell = True)
        subprocess.call('rm '+Workspace_laptop+'*',shell =True)
    except:
        print Root_laptop
        print 'error'
        sys.exit()
#On recuppere les fichier sur malbec
try:
    subprocess.call('rsync -zav -e"ssh -p11270" thorey@localhost:' +Root_malbec+ ' ' + Root_laptop[:-1], shell = True)
except:
    print 'error'
    sys.exit()

# On modifie le file E_Load.py pour le nom du dossier
with open( str(Routine_python)+'.py', 'r') as script:
    with open(str(Routine_python)+'_tmp.py', 'wr+') as script_tmp:
        for l in script:
            if l == 'c_input =  Name_Run_Folder\n':
                to_write = l.replace('Name_Run_Folder','"'+str(Root_laptop)+'"')
            elif l == 'c_output = Name_Workspace_Folder\n':
                to_write = l.replace('Name_Workspace_Folder','"'+str(Workspace_laptop)+'"')
            else:
                to_write = l
            script_tmp.write(to_write)
            
# On le lance
try:
    subprocess.call('python '+str(Routine_python)+'_tmp.py', shell = True)
except:
    print 'error'
    sys.exit()
