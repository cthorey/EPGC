import os
Root = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/Code_ELAS/'
Racine = Root
env = Environment()          
List_file = [elt for elt in os.listdir(Racine) if (elt.split('.')[-1] == 'f90') and (elt != 'Module_Init.f90') and (elt != 'Module_Init_tmp.f90')]        
env.Program(target = 'G1',source = List_file)

