import subprocess
import sys
import os
import shutil
import pydot
from struct import pack
import numpy as np
import pandas as pd
import seaborn as sns

def m_Max_Time(root,run,workspaces):
    input_file = root+'/'+workspaces[run]
    with open(input_file, 'rb') as f:
        try:
            data_pickle = pickle.load(f)
        except:
            data_pickle = pd.read_pickle(input_file)
    return '%3.1f'%(pd.DataFrame(data_pickle['Max']).tm.max())

def Bug_tcheque(run,Bugs):
    if run in Bugs:
        return True
    else:
        return False

def Map_Old_New_Name(run):

    def mapping_1(f):
        return f[0]+str(float(f[1:].replace('D','e')))
    def mapping_2(f):
        return f[0:2]+str(float(f[2:].replace('D','e')))
    def mapping_3(f):
        return f[0:1]+str(1./float(f[1:].replace('D','e')))    
    bug = run.split('_')[:-1]
    bug_n = []
    if len([f for f in bug if f[0]=='G']) == 1:
        bug_n.append(mapping_1(bug[0]).replace('E','el'))
        bug_n.append(mapping_1(bug[1]).replace('G','gr'))
        bug_n.append(mapping_1(bug[4]).replace('D','de'))
        bug_n.append(mapping_1(bug[7]).replace('S','si'))
        bug_n.append(mapping_1(bug[2]).replace('N','nu'))
        bug_n.append(mapping_3(bug[3]).replace('P','Pe'))
        bug_n.append(mapping_1(bug[5]).replace('C','Ps'))
        bug_n.append(mapping_1(bug[6]).replace('R','N1'))
        bug_n.append(mapping_2(bug[10]).replace('Dt','Dt'))
        bug_n.append(mapping_2(bug[8]))
        bug_n.append(mapping_2(bug[9]).replace('Ep','ep'))
    else:
        bug_n.append(mapping_1(bug[0]).replace('E','el'))
        bug_n.append(mapping_1(bug[1]).replace('G','gr'))
        bug_n.append(mapping_1(bug[4]).replace('D','de'))
        bug_n.append(mapping_1(bug[9]).replace('S','si'))
        bug_n.append(mapping_1(bug[2]).replace('N','nu'))
        bug_n.append(mapping_3(bug[3]).replace('P','Pe'))
        bug_n.append(mapping_1(bug[5]).replace('C','Ps'))
        bug_n.append(mapping_1(bug[6]).replace('R','N1'))
        bug_n.append(mapping_1(bug[7]).replace('G','ga'))
        bug_n.append(mapping_1(bug[8]).replace('I','In'))
        bug_n.append(mapping_2(bug[12]).replace('Dt','Dt'))
        bug_n.append(mapping_2(bug[10]))
        bug_n.append(mapping_2(bug[11]).replace('Ep','ep'))
    bug_new = '_'.join(bug_n)
    
    return bug_new

def Draw_Tree(Root,Startup_Node): 

    List = [f for f in os.listdir(Root) if f[:2] == 'el']
    Runs = ['_'.join([g for g in f.split('_') if g[0]!='M']) for f in List]  # On enleve M
    workspaces = dict(zip(Runs,List))
    Bugs = [f for f in os.listdir(Root+'/Run_Code') if f[0] == 'E' and f.endswith('BUG')]
    Bugs = map(Map_Old_New_Name,Bugs)
    graph = pydot.Dot(graph_type='graph')
    Level = [f[:2] for f in Runs[0].split('_')]
    color = ['#'+pack("BBB",*tuple(np.array(triplet)*255)).encode('hex') for triplet in sns.color_palette('deep',len(Level))]
    def Visit_Branch(List_Run,Parent_Node):
        for key_node,label_node in Parent_Node.iteritems():
            Level_Parent = Level.index(key_node[:2])
            Child_List = list(set([key.split('_')[Level_Parent+1] for key in List_Run if key.split('_')[Level_Parent] == key_node]))
            Tree_Child = [f for f in List_Run if f.split('_')[Level_Parent] == key_node]
            Child = dict(zip(Child_List,[label_node+'_'+f for f in Child_List]))
            node_parent = pydot.Node(label_node, label= key_node,style="filled", fillcolor=color[Level_Parent])
            graph.add_node(node_parent)
            for key_child,label_child in Child.iteritems():
                node_child = pydot.Node(label_child,label = key_child)
                graph.add_node(node_child)
                graph.add_edge(pydot.Edge(node_parent, node_child))
            if not Child_List[0][:2]  == 'ep':
                Visit_Branch(Tree_Child,Child)
            else:
                for key_child,label_child  in Child.iteritems():                
                    node_final_l = label_node+key_child
                    print Child.values()[0] in Bugs
                    if Child.values()[0] in Bugs:
                        color_Bug = 'red'
                    else:
                        color_Bug = 'white'
                    node_final = pydot.Node(node_final_l,
                                            label = str(m_Max_Time(Root,label_child,workspaces)),
                                            style="filled",
                                            fillcolor=color_Bug)
                    graph.add_node(node_final)
                    graph.add_edge(pydot.Edge(label_child, node_final))
                    print m_Max_Time(Root,label_child,workspaces)

    Visit_Branch(Runs,{'el0.0':'el0.0'})
    graph.write_pdf(Root+'/Tree_Runs.pdf')

Root = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/GRAV/Workspace_2015-04-21_0'
Startup_Node = {'el0.0':'el0.0'}
Draw_Tree(Root,Startup_Node)

