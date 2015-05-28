import os

Root = '/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/Code/EPGC/'
Racine = Root
env = Environment(CC = 'gfortran',CCFLAGS = '-O0')
# env.Append(CCFLAGS = ['-g','-O0','ftz','-ftree-vectorizer-verbose=5'])

List_file = [elt for elt in os.listdir(Racine) if (elt.split('.')[-1] == 'f90') and (elt != 'Module_Init.f90') and (elt != 'Module_Init_tmp.f90')]        
env.Program(target = 'run',source = List_file)


##### OPTION
##### -O3 -fp-model fast=2 ! Faster option
##### But first test -O0, the safer option, always use -g

