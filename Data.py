# Recupere les donne sur malbec,
# Load les donnes

import subprocess
import sys
import os

Folder_malbec = 'Run_New_4/'
Folder_laptop = 'Run_New_4/'
Workspace_laptop = 'Workspace_New_4/'
Root_malbec = '/gpfs/users/thorey/ELAS/'+Folder_malbec
Root_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/ELAS/'+Folder_laptop
Workspace_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/ELAS/'+Workspace_laptop
Routine_python = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Routine/ELAS/E_Load.py'

try :
    subprocess.call('rm -R '+Root_laptop +'*',shell = True)
    subprocess.call('rm '+Workspace_laptop+'*',shell =True)
except:
    print Root_laptop
    print 'error'
    sys.exit()
try:
    subprocess.call('rsync -zav -e"ssh -p11270" thorey@localhost:' +Root_malbec+ ' ' + Root_laptop, shell = True)
except:
    print 'error'
    sys.exit()

try:
    subprocess.call('python '+Routine_python, shell = True)
except:
    print 'error'
    sys.exit()
