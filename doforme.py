#!/usr/bin/python

import glob
import os
import sys

def writeVa(jobname):
    with open("va.pbs", 'w') as fid:
        fid.write("""####  PBS preamble

#PBS -N %s
#PBS -M chaomy@umich.edu
#PBS -m e

#PBS -l nodes=1:ppn=12:ib,pmem=4gb,walltime=240:00:00
#PBS -j oe
#PBS -V

#PBS -A qiliang_flux
#PBS -q flux
#PBS -l qos=flux

####  End PBS preamble

cd $PBS_O_WORKDIR

main()
{
    mpirun  vasp > Log
}

time main;
    """%(jobname))
        fid.close()
    return


def doforme():
    DirList  = ['FCC', 'Hcp', 'V0.9', 'V1.1', 'T1200', 'T2200', 'T5200', 'Vacancy']
    Root = os.getcwd()
    for Dir in DirList:
        os.chdir(Dir)

        if not os.path.isdir("New_continue"):
            os.mkdir("New_continue")

        os.system("cp Continue/* New_continue")
        os.system("cp CONTCAR  New_continue/POSCAR")
        os.system("cp ../INCAR New_continue")
        writeVa(Dir)
        os.system("mv va.pbs New_continue")
        os.chdir("New_continue")
     #   os.system("qsub va.pbs")


        os.chdir(Root)
    return

def collectConfig_md():
    DirList  = ['Perf', 'FCC', 'Hcp', 'V0.9', 'V1.1', 'T1200', 'T2200', 'T5200', 'Vacancy',
            'Monos_-0.01', 'Monos_0.01',  'Monos_-0.02', 'Monos_0.02',
            'Othos_-0.01', 'Othos_0.01',  'Othos_-0.02', 'Othos_0.02']

    Root = os.getcwd()
    for Dir in DirList:
        os.chdir(Dir)
        ConfigName = "Config-%s-3000"%(Dir)
        os.system("vasp2force -s 3000 > %s"%(ConfigName))
        os.system("cp  %s /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_3000"%(ConfigName))

        os.chdir(Root)
    return

def collectConfig_def(Num):
    DirList  = ['Perf', 'FCC', 'Hcp', 'V0.9', 'V1.1', 'T1200', 'T2200', 'T5200', 'Vacancy',
                'Monos_-0.01', 'Monos_0.01',  'Monos_-0.02', 'Monos_0.02',
                'Othos_-0.01', 'Othos_0.01',  'Othos_-0.02', 'Othos_0.02']
                #'Tpath_0.05',  'Tpath_0.15',  'Opath_0.05',  'Opath_0.15']

    Root = os.getcwd()
    os.system(": > /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_%d/dummy.config"%(Num))
    for Dir in DirList:
        os.chdir(Dir)
        os.chdir("Continue")
        ConfigName = "Config-%s-%d"%(Dir,Num)
        os.system("vasp2force -s %d > %s"%(Num, ConfigName))
        os.system("cp  %s /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_%d"%(ConfigName, Num))
        os.system("cat %s >> /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_%d/dummy.config"%(ConfigName, Num))
        os.chdir(Root)
    return

def collectConfig_final():
    DirList  = ['Perf', 'FCC', 'Hcp', 'V0.9', 'V1.1', 'T1200', 'T2200', 'T5200', 'Vacancy',
                'Monos_-0.01', 'Monos_0.01',  'Monos_-0.02', 'Monos_0.02',
                'Othos_-0.01', 'Othos_0.01',  'Othos_-0.02', 'Othos_0.02',
                'Tpath_0.05',  'Tpath_0.15',  'Opath_0.05',  'Opath_0.15']

    Root = os.getcwd()
    os.system(": > /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_final/dummy.config")
    for Dir in DirList:
        os.chdir(Dir)
        os.chdir("Continue")
        ConfigName = "Config-%s-final"%(Dir)
        os.system("vasp2force -f > %s"%(ConfigName))
        os.system("cp  %s /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_final"%(ConfigName))
        os.system("cat %s >> /scratch/qiliang_flux/chaomy/JOB/Nb/Fitapot/FITPOT_final/dummy.config"%(ConfigName))
        os.chdir(Root)
    return

def LoopIt():
    NumList = [3, 7, 9, 11, 13, 15, 17, 19, 20]
    for i in range(len(NumList)):
        Num = NumList[i]
        collectConfig_def(Num)
    return

if __name__ == '__main__':
  #  collectConfig_def(5)
  # LoopIt()
    collectConfig_final()

#cd New_Continue/
#cp ../Continue/* .
#cp ../CONTCAR ./POSCAR
