# Recupere les donne sur malbec,
# Load les donnes

import subprocess
import sys
import os

Folder_malbec = 'Run_New/'
Folder_laptop = 'Run_New/'
Root_malbec = '/gpfs/users/thorey/ELAS/'+Folder_malbec
Root_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/ELAS/'+Folder_laptop
Routine_python = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Routine/ELAS/E_Load.py'
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
