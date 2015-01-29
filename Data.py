#######################
# Recupere les donne sur malbec,
# Load les donnes

import subprocess
import sys
import os
import shutil
import pydot
import seaborn as sns
from struct import pack
import numpy as np
import pandas as pd

reset = 0
Folder = '2015-01-27_0'
REGIME = 'ELAS_GRAV'
Folder_malbec = 'Run_'+Folder+'/'
Folder_laptop = 'Run_'+Folder+'/'
Workspace_laptop = 'Workspace_'+Folder+'/'
Routine_python = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Routine/ELAS/E_Load'
Root_malbec = '/gpfs/users/thorey/'+REGIME+'/'+Folder_malbec
Root_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/'+REGIME+'/'+Folder_laptop
Workspace_laptop = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/'+REGIME+'/'+Workspace_laptop

#######################
# Fonction

def Tcheque_bug_file(source):
    print  'Les run qui ont bugger sont:',[ f for f in os.listdir(source+'Run_code/') if 'BUG' in f.split('_')]

def Tcheque_Number_file(source):
    Folder = [ f for f in os.listdir(source) if f[0]=='E']
    print 'On enleve :'
    print [f for f in Folder if len([g for g  in os.listdir(source+f) if g[0]=='R'])<3]
    map(shutil.rmtree,[source+f for f in Folder if len([g for g  in os.listdir(source+f) if g[0]=='R'])<3])

#######################
# main
    
# On fait un peu de menage
if reset == 1: # On enleve le dossier et tous ls workspace pr recommecner de zero
    try :
        subprocess.call('rm -R '+Root_laptop +'*',shell = True)
        subprocess.call('rm '+Workspace_laptop+'*',shell =True)
    except:
        print Root_laptop
        print 'error'
        sys.exit()
# On recuppere les fichier sur malbec
try:
    subprocess.call('rsync -zav -e"ssh -p11270" thorey@localhost:' +Root_malbec+ ' ' + Root_laptop[:-1], shell = True)
except:
    print 'error'
    sys.exit()

# On affiche le nombre de bug
Tcheque_bug_file(Root_laptop)

# On enleve les dossier qui vont faire bugger le truc
Tcheque_Number_file(Root_laptop)

#Tester si le dossier workspace existe
if not os.path.exists(Workspace_laptop):
    os.mkdir(Workspace_laptop)

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

# Un petit script pour afficher l'ensemble des run fait dans un graph !

def m_Max_Time(input_file):
        with open(input_file, 'rb') as f:
            try:
                data_pickle = pickle.load(f)
            except:
                data_pickle = pd.read_pickle(input_file)
        return pd.DataFrame(data_pickle['Max']).tm.max()

List = [f[:-4] for f in os.listdir(Workspace_laptop)]
List_Run = ['_'.join(f.split('_')[:8]) for f in List]
dict_Run = dict(zip(List_Run,os.listdir(Workspace_laptop)))
graph = pydot.Dot(graph_type='graph')
Level = ['E','G','N','P','D','C','R','S']
color = ['#'+pack("BBB",*tuple(np.array(triplet)*255)).encode('hex') for triplet in sns.color_palette('deep',8)]

def Visit_Branch(List_Run,Parent_Node):

    for key_node,label_node in Parent_Node.iteritems():
        Level_Parent = Level.index(key_node[0])
        Child_List = list(set([key.split('_')[Level_Parent+1] for key in List_Run if key.split('_')[Level_Parent] == key_node]))
        Tree_Child = [f for f in List_Run if f.split('_')[Level_Parent] == key_node]
        Child = dict(zip(Child_List,[label_node+'_'+f for f in Child_List]))
        node_parent = pydot.Node(label_node, label= key_node,style="filled", fillcolor=color[Level_Parent])
        graph.add_node(node_parent)
        for key_child,label_child in Child.iteritems():
            node_child = pydot.Node(label_child,label = key_child)
            graph.add_node(node_child)
            graph.add_edge(pydot.Edge(node_parent, node_child))
        if not Child_List[0][0]  == 'S':
            Visit_Branch(Tree_Child,Child)
        else:
            node_final_l = label_node+Child.keys()[0]
            node_final = pydot.Node(node_final_l,label = str(m_Max_Time(Workspace_laptop+'/'+dict_Run[Child.values()[0]])))
            graph.add_node(node_final)
            graph.add_edge(pydot.Edge(Child.values()[0], node_final))
            print m_Max_Time(Workspace_laptop+'/'+dict_Run[Child.values()[0]])
if REGIME == 'GRAV':
    Visit_Branch(List_Run,{'E0.0':'E0.0'})
else:
    Visit_Branch(List_Run,{'E1.0':'E1.0'})
graph.write_pdf(Workspace_laptop+'/GRAPH_'+REGIME+'.pdf')
