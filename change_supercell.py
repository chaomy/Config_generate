#!/usr/local/bin/python
# encoding: utf-8
#
######################################################
######## This is transform the Bcc 2 atom ############
######## to Bcc one atom config program   ############
######################################################

import  os, shutil
import  glob, re
import  numpy as np

class  ChangeSuperCell(object):
    def __init__(self, Num):
        self._Element = 'WRe',
        self._Num = Num
        self._LatticeList = [5.98753,
                             5.98281,
                             5.97809,
                             5.97337,
                             5.96924,
                             5.94977]
        self._Kpoints = [33, 33, 33]
        self._effectiveMass = [(183.84 * (1 - 0.05 * i) + 186.21 *
            0.05 * i) for i in range(1,6)]
        self._effectiveMass.append((183.84 * 0.5 + 186.21 * 0.5))
        self._Potentialist = ['WRe.0-05.fhi.UPF',
        'WRe.0-1.fhi.UPF',
        'WRe.0-15.fhi.UPF',
        'WRe.0-2.fhi.UPF',
        'WRe.0-25.fhi.UPF',
        'WRe.0-5.fhi.UPF']
        return

    def gnInfile(self, CellMatrix, filename):
        with open(filename,'w') as fid:
            fid.write("""
&control
calculation='scf',
prefix= 'WRe',
outdir='./results',
tstress = .true.,
tprnfor = .true.,
pseudo_dir = './',
etot_conv_thr=1.0D-5,
forc_conv_thr=1.0D-4,
/
&system
    ibrav= 0, nat=  1, ntyp= 1,
    occupations='smearing',
    smearing='m-p',
    degauss=0.02D0,
    ecutwfc =45.0,
/
&electrons
conv_thr    = 1.D-10,
/
&ions
ion_dynamics='bfgs',
/
&cell
cell_dynamics='bfgs',
press_conv_thr=0.1D-0,
/
CELL_PARAMETERS {bohr}
%f  %f  %f
%f  %f  %f
%f  %f  %f
ATOMIC_SPECIES
%s  %f  %s
ATOMIC_POSITIONS {bohr}
%s   0  0  0
K_POINTS automatic
%d %d %d  0 0 0
"""%( CellMatrix[0,0],
        CellMatrix[0,1],
        CellMatrix[0,2],

        CellMatrix[1,0],
        CellMatrix[1,1],
        CellMatrix[1,2],

        CellMatrix[2,0],
        CellMatrix[2,1],
        CellMatrix[2,2],

        'W',
        self._effectiveMass[self._Num],
        self._Potentialist[self._Num],

        'W',
        self._Kpoints[0],
        self._Kpoints[1],
        self._Kpoints[2],
        ))
        return

    def getMatrix(self,
                  infile):
        with open(infile, 'r') as fid:
            Raw = fid.read()
            fid.close()

        realN = r'(-?\d*\.\d*)'
        In = r'\s*'
        lineHead = r'CELL_PARAMETERS\s*\{bohr\}\s*'
        lineNumber = realN + In + \
                    realN + In + \
                    realN + In
        M = lineHead + \
                lineNumber + lineNumber + lineNumber
        findPosition = re.compile(M,
                                re.DOTALL)
        PositionM = findPosition.findall(Raw)
        if PositionM:
            print PositionM

        CellMatrix = np.zeros([3,3],"float")
        for i in range(3):
            for j in range(3):
                CellMatrix[i,j] = float(PositionM[0][i + j * 3])
        return np.mat(CellMatrix)

    def MyMatrixMul(self, m1, m2):
        m = np.zeros([3,3],'float')
        for i in range(3):
            m[i, 0] = m1[i, 0]* m2[0, 0] + m1[i, 1] * m2[1, 0] + m1[i, 2]*\
                    m2[2, 0]
            m[i, 1] = m1[i, 0]* m2[0, 1] + m1[i, 1] * m2[1, 1] + m1[i, 2]*\
                    m2[2, 1]
            m[i, 2] = m1[i, 0]* m2[0, 2] + m1[i, 1] * m2[1, 2] + m1[i, 2]*\
                    m2[2, 2]
        m = np.mat(m)
        return m

    def QE_Bcc_two_to_One(self, filename):
        Newfilename = 'new-%s'%(filename)

        Lattice = self._LatticeList[5 - 1]

        OldMatrix = self.getMatrix(filename)
        OldStrain = OldMatrix/ Lattice

        New_Cordi_M = np.mat([[-0.5, 0.5, 0.5],
                              [0.5, -0.5, 0.5],
                              [0.5, 0.5, -0.5]], 'float')

        #NewStrain = OldtoNew_M * OldStrain * OldtoNew_M.transpose()

        NewStrain= OldStrain

       # Sv1 = self.MyMatrixMul(OldtoNew_M.transpose(), OldStrain)
       # Sv2 = self.MyMatrixMul(Sv1, OldtoNew_M)
       # print "Use my", Sv2

        NewMatrix = NewStrain * New_Cordi_M.transpose()
        NewMatrix = NewMatrix.transpose()

        print NewMatrix
        NewMatrix = NewMatrix * Lattice
        self.gnInfile(NewMatrix , Newfilename)
        return

    def LoopChangeInfiles(self):
        FileList = glob.glob("WRe.in.*")
        for file in FileList:
            self.QE_Bcc_two_to_One(file)

        return

if __name__ == '__main__':
    A = ChangeSuperCell(4)
    A.LoopChangeInfiles()








