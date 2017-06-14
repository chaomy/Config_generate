#!/usr/bin/python 

import glob 
import os 
import sys 

def Combine_configs(): 
    DirList  = ['Perf', 'FCC', 'Hcp', 'V0.9', 'V1.1', 'T1200', 'T2200', 'T5200', 'Vacancy',
            'Monos_-0.01', 'Monos_0.01',  'Monos_-0.02', 'Monos_0.02',
            'Othos_-0.01', 'Othos_0.01',  'Othos_-0.02', 'Othos_0.02'] 
    Step = '10'  

    ConfigList = [] 
    for i in range(DirList): 
        ConfigList.append('%s-%s'%(DirList[i],Step)) 

    os.system(": > dummy.config") 
    List = [1,2,3,4,5]
    for i in range(List):
        os.system("cat %s >> dummy.config"%(ConfigList[List[i]]) 

if __name__ == '__main__': 
    Combine_configs() 
