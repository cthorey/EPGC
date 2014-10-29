import subprocess
import sys
import os

try:
    subprocess.call('rm Test/Run/*', shell = True)
except:
    print 'error'
    sys.exit()
try:
    subprocess.call('rm Test/Workspace/*', shell = True)
except:
    sys.exit()
try:
    subprocess.call('scons', shell = True)
except:
    sys.exit()
try:
    subprocess.call('./G1', shell = True)
except:
    pass
try:
    subprocess.call('python ../Routine/TEST/Load.py', shell=True)
except:
    sys.exit()
