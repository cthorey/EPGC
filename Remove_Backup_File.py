import os,sys

Racine = '//gpfs/users/thorey'
folders = []
folders.append(os.path.join(Racine,'ELAS','MSkin_TSc_Newton_HSc_Newton_RBercovici','Run_2015-05-30_0'))
folders.append(os.path.join(Racine,'ELAS','MSkin_TSc_Newton_HSc_Newton_RBercovici','Run_2015-06-02_0'))
folders.append(os.path.join(Racine,'ELAS','MSkin_TSc_Newton_HSc_Newton_RArrhenius','Run_2015-06-02_0'))
folders.append(os.path.join(Racine,'ELAS','MSkin_TSc_Newton_HSc_Newton_RArrhenius','Run_2015-06-11_0'))
folders.append(os.path.join(Racine,'GRAV','MSkin_TSc_GFD_HSc_Newton_RArrhenius','Run_2015-05-30_0'))
folders.append(os.path.join(Racine,'GRAV','MSkin_TSc_GFD_HSc_Newton_RBercovici','Run_2015-06-23_3'))
folders.append(os.path.join(Racine,'GRAV','MSkin_TSc_GFD_HSc_Newton_RArrhenius','Run_2015-06-23_2'))
folders.append(os.path.join(Racine,'GRAV','MSkin_TSc_GFD_HSc_Newton_RArrhenius','Run_2015-06-24_0'))
folders.append(os.path.join(Racine,'GRAV','MSkin_TSc_GFD_HSc_Newton_RArrhenius','Run_2015-06-24_1'))
folders.append(os.path.join(Racine,'GRAV','MSkin_TSc_GFD_HSc_Newton_RBercovici','Run_2015-06-24_0'))
folders.append(os.path.join(Racine,'ELASGRAV','MSkin_TSc_Newton_HSc_Newton_RBercovici','Run_2015-07-05_0'))
folders.append(os.path.join(Racine,'ELASGRAV','MSkin_TSc_Newton_HSc_Newton_RArrhenius','Run_2015-06-26_0'))
def Extract_Compteur(s):
    return int(s.split('_')[1].split('.')[0])

for folder in folders:
    print folder
    for i,f in enumerate(os.walk(folder)):
        if len(f[2]) != 0:
            Backup_file = [g for g in os.listdir(f[0]) if g.split('_')[0] == 'Backup']
            try:
                Last_Rv_file = Extract_Compteur([g for g in os.listdir(f[0]) if g.split('_')[0] == 'RV'][-1])
                Rm_Backup = [os.path.join(f[0],g) for g in Backup_file if Extract_Compteur(g)<Last_Rv_file-10]
                map(os.remove,Rm_Backup)
            except:
                print 'Ce dossier a buger %s'%(str(f[0]))

            
