#!/usr/bin/python  

import glob 
import os 

def doforme():
    DirList  = ['FCC', 'Hcp', 'V0.9', 'V1.1', 'T1200', 'T2200', 'T5200', 'Vacancy'] 
    DirList_1 = []

    Root = os.getcwd() 
    for Dir in DirList:  
        FileName = "Out_%s"%(DirList) 
        os.chdir(Dir) 

        if os.path.isdir("New_Continue"):
            os.system("cp New_Continue/OUTCAR  ../%s"%(FileName)) 

        os.chdir(Root) 
    return 

if __name__ == '__main__':
    doforme() 
    

#cd New_Continue/
#cp ../Continue/* .
#cp ../CONTCAR ./POSCAR  
