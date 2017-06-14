#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./gn_config.py
#
###################################################################
#
# Purpose :  generate configurations for LMP and VASP
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################

import os
import numpy as np
import ase.lattice.cubic as Cubic
import ase.lattice.hexagonal as Hexagonal
import ase.io
from ase.lattice.spacegroup import crystal
import copy


class gnStructure(object):

    def __init__(self, pot=None):
        self._cartition = [[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1]]

        self._default_direction = [[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]

        print pot
        self._default_size = (5, 5, 5)
        self._config_file_format = 'vasp'
        self._lattice_constant = pot['lattice']
        self._element = pot['element']
        return

    def mymkdir(self, dirname):
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        return

    def write_config_with_fix(self, atoms):
        #  a_list = np.array([102, 107, 103, 105, 104]);
        #  b_list = np.array([123, 122, 124, 84,   89]);
        fix_list = np.array([107, 103, 105, 122, 124, 84])

        cell = atoms.get_cell()
        positions = atoms.get_positions()
        print positions
        with open("POSCAR", 'w') as fid:
            fid.write("W\n")
            fid.write("1\n")
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (cell[0, 0], cell[0, 1], cell[0, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (cell[1, 0], cell[1, 1], cell[1, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (cell[2, 0], cell[2, 1], cell[2, 2]))
            fid.write("W\n")
            fid.write("%d\n" % (len(atoms)))
            fid.write("Selective dynamics\n")
            fid.write("Cartesian\n")
            for i in range(len(atoms)):
                if ((i == fix_list).any()):
                    fid.write("%12.6f %12.6f %12.6f   T   T   F\n"
                              % (positions[i, 0], positions[i, 1], positions[i, 2]))
                else:
                    fid.write("%12.6f %12.6f %12.6f   T   T   T\n"
                              % (positions[i, 0], positions[i, 1], positions[i, 2]))
            fid.close()
        os.system("cp POSCAR POSCAR.vasp")
        return

    def vasp_intro_alloy(self):
        atoms = ase.io.read(filename="POSCAR",
                            format='vasp')
        print atoms
        alloy_num = 0
        for atom in atoms:
            if np.random.rand() < 0.1:
                atom.symbol = 'Re'
                alloy_num += 1

        cell = atoms.get_cell()
        with open("POSCAR", 'w') as fid:
            fid.write("W   Re\n")
            fid.write("1\n")
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (cell[0, 0], cell[0, 1], cell[0, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (cell[1, 0], cell[1, 1], cell[1, 2]))
            fid.write("%12.6f %12.6f %12.6f\n" %
                      (cell[2, 0], cell[2, 1], cell[2, 2]))
            fid.write("W   Re\n")
            fid.write("%d   %d\n" % ((len(atoms) - alloy_num), alloy_num))
            fid.write("Cartesian\n")
            for atom in atoms:
                if atom.symbol == 'W':
                    fid.write("%12.6f %12.6f %12.6f\n"
                              % (atom.position[0], atom.position[1], atom.position[2]))
            for atom in atoms:
                if atom.symbol == 'Re':
                    fid.write("%12.6f %12.6f %12.6f\n"
                              % (atom.position[0], atom.position[1], atom.position[2]))
            fid.close()
        os.system("cp POSCAR POSCAR.vasp")
        return

    def cut_plane_for_vacancy(self, atoms):
        cell = atoms.get_cell()
        print cell[2, 2]
        keep_atom = cell[2, 2] - 14.0
        index_list = []

        for i in range(len(atoms)):
            atom = atoms[i]
            if atom.position[2] > keep_atom:
                index_list.append(atom.index)
        del atoms[index_list]
        return atoms

    def lmp_change_box(self, in_cell):
        unit_x = np.linalg.norm(in_cell[0, :])
        unit_y = np.linalg.norm(in_cell[1, :])
        unit_z = np.linalg.norm(in_cell[2, :])

        vect_x = in_cell[0, :] / unit_x
        vect_y = in_cell[1, :] / unit_y
        vect_z = in_cell[2, :] / unit_z

        ax = unit_x
        bx = np.dot(in_cell[1, :], vect_x.transpose())
        by = np.linalg.norm(np.cross(vect_x, in_cell[1, :]))

        cx = np.dot(in_cell[2, :],  vect_x.transpose())

        a_cross_b = np.cross(in_cell[0, :], in_cell[1, :])
        a_cross_b_unit = a_cross_b / np.linalg.norm(a_cross_b)

        Acvect_ycvect_x = np.cross(a_cross_b_unit, vect_x)

        cy = np.dot(in_cell[2, :], Acvect_ycvect_x.transpose())
        cz = np.abs(np.dot(in_cell[2, :], a_cross_b_unit.transpose()))

        out_cell = np.mat(np.zeros([3, 3], "float"))

        out_cell[0, 0], out_cell[0, 1], out_cell[0, 2] = ax, 0,  0
        out_cell[1, 0], out_cell[1, 1], out_cell[1, 2] = bx, by, 0
        out_cell[2, 0], out_cell[2, 1], out_cell[2, 2] = cx, cy, cz
        return out_cell

    #################################################################
    # Sun Mar 26 13:00:28 2017
    #################################################################
    def write_lmp_config_data(self,
                              atoms,
                              filename="lmp_init.txt"):

        positions = atoms.get_positions()
        atom_num = len(positions)
        cell = atoms.get_cell()
        print "cell is "
        print cell

        with open(filename, mode="w") as fout:
            fout.write("#lmp data config")
            fout.write("\n")
            fout.write("%d atoms\n" % (atom_num))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n" % (0, cell[0, 0]))
            fout.write("%f\t%f ylo yhi\n" % (0, cell[1, 1]))
            fout.write("%f\t%f zlo zhi\n" % (0, cell[2, 2]))
            fout.write("%f  %f  %f xy xz yz\n"
                       % (cell[1, 0],
                          cell[2, 0],
                          cell[2, 1]))

            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(atom_num):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"
                           % (i + 1, positions[i, 0], positions[i, 1], positions[i, 2]))
        fout.close()
        return

    def write_lmp_coords(self, atoms, filename="lmp_coord"):
        positions = atoms.get_positions()
        atom_num = len(positions)

        with open(filename, mode="w") as fid:
            fid.write("%d\n" % (atom_num))
            for i in range(atom_num):
                fid.write("%d  %f  %f  %f \n"
                          % (i + 1, positions[i, 0], positions[i, 1], positions[i, 2]))
            fid.close()
        return

    def write_poscar(self, atoms, filename="POSCAR"):
        ase.io.write(filename=filename,
                     images=atoms,
                     format='vasp')
        os.system("cp POSCAR POSCAR.vasp")
        return

    def write_lmp_cfg(self, atoms,
                      filename=None):
        if filename is None:
            locfilename = "lmp_init.cfg"
        else:
            locfilename = filename
        ase.io.write(filename=locfilename,
                     images=atoms,
                     format='cfg')
        return

    def set_config_file_format(self, in_format):
        self._config_file_format = in_format
        return

    def write_config_output(self, atoms):
        if self._config_file_format == 'vasp':
            self.write_poscar(atoms)
        elif self._config_file_format == 'lmp':
            self.write_lmp_cfg(atoms)
        return

    def set_directions(self, direction):
        self._directions = direction
        return

    def set_element(self, symbol):
        self._element = symbol
        return

    def set_lattce_constant(self, lattice_constant):
        self._lattice_constant = lattice_constant
        return

    def set_size(self, in_size):
        self._default_size = in_size
        return


class add_strain(object):

    def __init__(self):
        self._delta = 0.02
        return

    def volumetric_strain(self, delta=None):
        ###### c11 #####
        if delta == None:
            delta = self._delta
        strain = np.mat([[1 + delta, 0, 0],
                         [0, 1 + delta, 0],
                         [0, 0, 1 + delta]], "float")
        return strain

    def volume_conserving_orthohombic(self, delta=None):
        ###### c12 #####
        if delta is None:
            delta = self._delta
        strain = np.mat([[1 + delta, 0, 0],
                         [0, 1 - delta, 0],
                         [0, 0, 1 / (1 - delta**2)]], 'float')
        return strain

    def volume_conserving_monoclinic(self, delta=None):
        ###### c44 #####
        if delta is None:
            delta = self._delta
        strain = np.mat([[1, delta, 0],
                         [delta, 1, 0],
                         [0, 0, 1 / (1 - delta**2)]], 'float')
        return strain


class bcc(gnStructure, add_strain):

    def __init__(self, pot=None):
        self.pot = pot
        print self.pot
        gnStructure.__init__(self, self.pot)
        add_strain.__init__(self)
        self._lattice_constant = self.pot['latbcc']
        self._element = pot['element']
        return

    def set_bcc_primitive_direction(self):
        self._primitive_directions = np.array([[-0.5, 0.5, 0.5],
                                               [0.5, -0.5, 0.5],
                                               [0.5,  0.5, -0.5]])
        return

    def set_bcc_primitive(self, in_size=(1, 1, 1)):
        atoms = ase.atoms.Atoms(symbols=self._element,
                                positions=[(0, 0, 0)],
                                info={'unit_cell': 'primitive'},
                                pbc=(1, 1, 1))

        self.set_bcc_primitive_direction()
        cell = self._primitive_directions * self._lattice_constant
        atoms.set_cell(cell, scale_atoms=True)
        atoms = atoms.repeat(in_size)
        return atoms

    def write_bcc_primitive(self, in_size=None):
        atoms = self.set_bcc_primitive(in_size)
        self.write_config_output(atoms)
        return

    def set_bcc_convention(self,
                           in_direction=None,
                           in_size=(1, 1, 1)):
        if in_direction is None:
            in_direction = self._default_direction
        atoms = Cubic.BodyCenteredCubic(directions=in_direction,
                                        latticeconstant=self._lattice_constant,
                                        size=in_size,
                                        symbol=self._element,
                                        pbc=(1, 1, 1))
        return atoms

    def write_bcc_convention(self, in_direction, in_size=(1, 1, 1)):
        atoms = self.set_bcc_convention(in_direction, in_size)
        self.write_config_output(atoms)
        return

    def write_bcc_with_strain(self,
                              delta=None,
                              in_tag='ortho',
                              in_size=None):

        if in_tag == 'ortho' or in_tag == 'c12':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono' or in_tag == 'c44':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume' or in_tag == 'c11':
            strain = self.volumetric_strain(delta)

        atoms = self.set_bcc_convention(in_size)
        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)
        self.write_config_output(atoms)
        return

    def write_bcc_primitive_with_strain(self,
                                        delta=None,
                                        in_tag=None,
                                        in_size=None):
        if in_tag == 'ortho' or in_tag == 'c12':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono' or in_tag == 'c44':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume' or in_tag == 'c11':
            strain = self.volumetric_strain(delta)

        atoms = self.set_bcc_primitive(in_size)

        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)

        self.write_config_output(atoms)
        return

    def example_of_crystal(self):
        a = 3.3
        Nb = crystal('Nb',
                     [(0, 0, 0)],
                     spacegroup=229,
                     cellpar=[a, a, a, 90, 90, 90],
                     primitive_cell=True)
        print Nb.get_cell()
        Nb.write('geometry_Bcc.in')
        return


class fcc(gnStructure, add_strain):

    def __init__(self, pot=None):
        gnStructure.__init__(self, pot)
        add_strain.__init__(self)
        if pot is not None:
            self._lattice_constant = pot['latfcc']
        return

    def set_fcc_primitive_direction(self):
        self._primitive_directions = np.array([[0.0, 0.5, 0.5],
                                               [0.5, 0, 0.5],
                                               [0.5, 0.5, 0.0]], 'float')
        return

    def set_fcc_primitive(self, in_size=None):
        if in_size == None:
            in_size = (1, 1, 1)
        atoms = ase.atoms.Atoms(symbols=self._element,
                                positions=[(0, 0, 0)],
                                info={'unit_cell': 'primitive'},
                                pbc=(1, 1, 1))

        self.set_fcc_primitive_direction()
        cell = self._primitive_directions * self._lattice_constant
        atoms.set_cell(cell, scale_atoms=True)
        atoms2 = atoms.repeat(in_size)
        return atoms2

    def set_fcc_convention(self,
                           in_direction=None,
                           in_size=(4, 4, 4)):
        l_size = in_size
        if in_direction is not None:
            l_direction = in_direction
        else:
            l_direction = self._default_direction

        atoms = Cubic.FaceCenteredCubic(directions=l_direction,
                                        latticeconstant=self._lattice_constant,
                                        size=l_size,
                                        symbol=self._element,
                                        pbc=(1, 1, 1))
        return atoms

    def write_fcc_convention(self, in_size=None):
        atoms = self.set_fcc_convention(in_size)
        self.write_config_output(atoms)
        return

    def write_fcc_primitive(self, in_size=None):
        atoms = self.set_fcc_primitive(in_size)
        self.write_config_output(atoms)
        return

    def write_fcc_with_strain(self,
                              delta=None,
                              in_tag='ortho',
                              in_size=None):
        if in_tag == 'ortho':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume':
            strain = self.volumetric_strain(delta)

        atoms = self.set_fcc_convention(in_size)

        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)

        self.write_config_output(atoms)
        return

    def write_fcc_primitive_with_strain(self,
                                        delta=None,
                                        in_tag=None,
                                        in_size=None):
        if in_tag == 'ortho' or in_tag == 'c12':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono' or in_tag == 'c44':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume' or in_tag == 'c11':
            strain = self.volumetric_strain(delta)

        atoms = self.set_fcc_primitive(in_size)

        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)

        self.write_config_output(atoms)
        return


class hcp(gnStructure, add_strain):

    def __init__(self, pot=None):
        gnStructure.__init__(self, pot)
        add_strain.__init__(self)
        if pot is not None:
            self.lattice_a = pot['ahcp']
            self.lattice_c = pot['chcp']
        return

    def set_hcp_lattice_constant_ratio(self,
                                       ta=np.sqrt(3) / 2.,
                                       tb=np.sqrt(8. / 3.)):
        self.lattice_a = self._lattice_constant * ta
        self.lattice_c = self.lattice_a * tb
        return

    def set_hcp_lattice_constant(self, a, c):
        self.lattice_a = a
        self.lattice_c = c
        return

    def set_hcp_direction(self):
        self.hcp_directions = np.array([[1, 0, 0],
                                        [-0.5, np.sqrt(3) / 2, 0],
                                        [0, 0, np.sqrt(8. / 3)]])
        return

    def set_hcp_convention(self, in_size=None):
        if in_size is None:
            in_size = (1, 1, 1)
        atoms = \
            Hexagonal.HexagonalClosedPacked(latticeconstant={'a': self.lattice_a,
                                                             'c': self.lattice_c},
                                            size=in_size,
                                            symbol=self._element,
                                            pbc=(1, 1, 1))
        print atoms.get_cell()
        return atoms

    def write_hcp_convention(self, size=None):
        atoms = self.set_hcp_convention(in_size=size)
        self.write_config_output(atoms)
        return

    def write_hcp_poscar(self, alat):
        with open("POSCAR", 'w') as fid:
            fid.write("""hcp
    %f
   1.00000         0.00000000000000      0.00000
  -0.50000         0.86602540378444      0.00000
   0.00000         0.00000000000000      1.7320
2
direct
   0.00000000000000    0.00000000000000      0.000000
   0.33333333333333    0.66666666666667      0.500000
            """ % (alat))
        return
