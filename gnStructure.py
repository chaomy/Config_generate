#!/usr/bin/python
# encoding: utf-8

import ase.lattice.cubic as Cubic
import copy
#from   ase.visualize import view
import ase.io
import numpy as np
import os
import re


class cal(object):

    def __init__(self,
                 directions,
                 size,
                 latticeconstant,
                 element,
                 structure):

        self._element = element
        self._directions = directions
        self._lattice_constant = latticeconstant
        self._structure = structure
        self._size = size

        self.exe = "mpirun vasp"

        self._atoms = self.get_SuperCell()
        self._PerfPosition = self._atoms.get_positions()
        self._PerfCells = self._atoms.get_cell()

    def get_SuperCell(self):
        if self._structure == 'bcc':
            atoms = Cubic.BodyCenteredCubic(directions=self._directions,
                                            size=self._size,
                                            latticeconstant=self._lattice_constant,
                                            symbol=self._element,
                                            pbc=(1, 1, 1))
        elif self._structure == 'fcc':
            atoms = Cubic.FaceCenteredCubic(directions=self._directions,
                                            size=self._size,
                                            latticeconstant=self._lattice_constant,
                                            symbol=self._element,
                                            pbc=(1, 1, 1))

        ase.io.write(filename="POSCAR",
                     images=atoms,
                     format='vasp')
        return atoms

    def gndisplacement(self,
                       Cell,
                       Positions,
                       DisplacementVector):
        AtomN = len(Positions)
        Displacement = copy.deepcopy(Positions)
        DisplacementVector = [0.4, 0.4, 0.0]

        for i in range(AtomN):
            if i < 0.5 * AtomN:
                for j in range(3):
                    Displacement[i][j] = 0.0
            elif i >= 0.5 * AtomN:
                Displacement[i] = DisplacementVector
        return Displacement

    def LoopGSF(self):
        for i in range(0, 20):
            DisplacementVector = [i * 0.05, 0, 0]
            DisplacementMatrix = self.gndisplacement(self._PerfCells,
                                                     self._PerfPosition,
                                                     DisplacementVector)

            Localatoms = self._atoms.copy()
            Localatoms.translate(DisplacementMatrix)

            ase.io.write(filename="POSCAR",
                         images=Localatoms,
                         format='vasp')
            os.system("%s" % (self.exe))
            (Energy, Stress, Volume) = self.VAgetData()
            self.output(0.05 * i,
                        Energy,
                        Stress)
            os.system("mv POSCAR POSCAR-{0:03d}.vasp".format(i))
        return

    def testGSF(self):
        for i in range(0, 20):
            DisplacementVector = [i * 0.05, 0, 0]
            DisplacementMatrix = self.gndisplacement(self._PerfCells,
                                                     self._PerfPosition,
                                                     DisplacementVector)

            Localatoms = self._atoms.copy()
            Localatoms.translate(DisplacementMatrix)

            ase.io.write(filename="POSCAR-{0:03d}.vasp".format(i),
                         images=Localatoms,
                         format='vasp')

        return

    def VAgetData(self):
        Stress = np.arange(6, dtype="float")
        Stress.shape = ([6, 1])
        with open("OUTCAR", 'r') as fid:
            real_N = r"(-?\d*\.\d{5})"
            real_Num = r"(-?\d*\.\d*)"
            In = r"\s*"
            for line in fid:
                matchV = \
                    re.search(r"\s*volume of cell : \s*" + real_Num,
                              line)
                matchE = \
                    re.search(r"\s*energy\s*without\s*entropy=\s*" +
                              real_Num + In + r"energy\(sigma->0\)\s*=\s*" +
                              real_Num, line)
                matchS = re.search(r"^\s*in kB\s*" + real_N +
                                   In + real_N + In + real_N + In + real_N +
                                   In + real_N + In + real_N, line)
                if matchS:
                    for i in range(6):
                        Stress[i] = float(matchS.group(i + 1))
                if matchE:
                    Energy = float(matchE.group(2))
                if matchV:
                    Volume = float(matchV.group(1))
            fid.close()
        return (Energy, Stress, Volume)

    def output(self,
               Displacement,
               Energy,
               Stress):
        with open("DATA", 'a') as fid:
            fid.write("%5.4f\t%10.8f\t%10.8f\t%10.8f\t"
                      "%10.8f\t%10.8f\t%10.8f\t%10.8f\n"
                      % (Displacement, Energy, Stress[0], Stress[1],
                         Stress[2], Stress[3], Stress[4], Stress[5]))
            fid.close()
        return


if __name__ == '__main__':
    Job = cal(directions=[[1, 0, 0],
                          [0, 1, 0],
                          [0, 0, 1]],
              size=(3, 3, 3),
              latticeconstant=4.2,
              element='Nb',
              structure='fcc')
    Job.get_SuperCell()
