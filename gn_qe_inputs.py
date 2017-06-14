#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name :
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

import numpy as np


class gn_qe_infile(object):

    def __init__(self,
                 pot):
        self.pot = pot
        self.cal_type = 'scf'
        self.degauss = None
        # set defaults
        self.set_degauss()
        self.set_cal_prefix()
        self.set_cal_type()
        self.set_ecut()
        self.set_disk_io()
        self.set_thr()
        return

    def set_cal_type(self, caltype='scf'):
        self.cal_type = caltype
        return

    def set_disk_io(self, diskinfo='none'):
        self.diskinfo = diskinfo
        return

    def set_degauss(self, degauss='0.2D0'):
        self.degauss = degauss
        return

    def set_cal_prefix(self, prefix='qe'):
        self.cal_prefix = prefix
        return

    def set_ecut(self, ecut='38'):
        self.ecutwfc = ecut
        return

    def set_thr(self, thr='1.0D-5'):
        self.conv_thr = thr
        return

    def qe_write_control(self, fid, atoms):
        fid.write("""&control
calculation = {},
prefix = {},
tstress = T
tprnfor = T
outdir = 'results',
pseudo_dir = './'
disk_io = {}
/
""".format(self.cal_type,
           self.cal_prefix,
           self.diskinfo))
        return fid

    def qe_write_system(self, fid, atoms):
        fid.write("""&system
ibrav= 0, nat= {}, ntyp = {},
occupations = 'smearing',
smearing= 'mp',
degauss = {},
ecutwfc = {},
/
""".format(atoms.get_number_of_atoms(),
           len(np.unique(atoms.get_chemical_symbols())),
            self.degauss,
            self.ecutwfc))
        return fid

    def qe_write_electrons(self, fid):
        fid.write("""&electrons
conv_thr = {},
/
&ions
ion_dynamics='bfgs',
/
""".format(self.conv_thr))
        return fid

    def gn_infile_dipole_ideal_shear(self,
                                     atoms=None):
        with open('qe.in', 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, (21, 21, 21))
            fid.close()
        return

    def gn_infile_dipole_screw_atoms(self,
                                     atoms=None):
        self.set_cal_type('relax')
        self.set_ecut('40')
        with open('{}.in'.format(self.pot['element']), 'w') as fid:
            fid = self.qe_write_control(fid, atoms)
            fid = self.qe_write_system(fid, atoms)
            fid = self.qe_write_electrons(fid)
            fid = self.qe_write_cell(fid, atoms.get_cell())
            fid = self.qe_write_species(fid, atoms, self.pot)
            fid = self.qe_write_pos(fid, atoms)
            fid = self.qe_write_kpts(fid, (1, 2, 8))
            fid.close()
        return

    def qe_write_cell(self, fid, cell):
        fid.write("CELL_PARAMETERS {angstrom}\n")
        fid.write("""{}  {}  {}
{}  {}  {}
{}  {}  {}
""".format(cell[0, 0], cell[0, 1], cell[0, 2],
           cell[1, 0], cell[1, 1], cell[1, 2],
           cell[2, 0], cell[2, 1], cell[2, 2]))
        return fid

    def qe_write_species(self, fid, atoms, pot):
        fid.write("ATOMIC_SPECIES\n")
        fid.write("{}  {}  {}\n".format(pot['element'],
                                        pot['mass'],
                                        pot['file']))
        return fid

    def qe_write_pos(self, fid, atoms):
        fid.write("ATOMIC_POSITIONS {alat}\n")
        sym = atoms.get_chemical_symbols()
        pos = atoms.get_scaled_positions()
        for i in range(len(pos)):
            fid.write("{}  {}  {}  {}\n".format(sym[i],
                                                pos[i, 0],
                                                pos[i, 1],
                                                pos[i, 2]))
        return fid

    def qe_write_kpts(self, fid, kpts=(1, 1, 1)):
        fid.write("K_POINTS automatic\n")
        fid.write("{:d} {:d} {:d}  0  0  0\n".format(
            kpts[0], kpts[1], kpts[2]))
        return fid
