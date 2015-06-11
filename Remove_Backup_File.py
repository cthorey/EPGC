import os,sys

Racine = '//gpfs/users/thorey'

folders.append()
folders = ['/Users/thorey/Documents/These/Projet/Refroidissement/Skin_Model/SCAPAD/ELAS/Run_New_8']

def Extract_Compteur(s):
    return int(s.split('_')[1].split('.')[0])

for i,f in enumerate(os.walk(folders[0])):
    if len(f[2]) != 0:
        Backup_file = [g for g in os.listdir(f[0]) if g.split('_')[0] == 'Backup']
        Last_Rv_file = Extract_Compteur([g for g in os.listdir(f[0]) if g.split('_')[0] == 'RV'][-1])
        Rm_Backup = [os.path.join(f[0],g) for g in Backup_file if Extract_Compteur(g)<Last_Rv_file-10]
        map(os.remove,Rm_Backup)

            
