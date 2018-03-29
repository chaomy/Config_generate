#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-28 21:31:42


import os
import numpy as np
import ase.lattice.cubic as Cubic
import ase.lattice.hexagonal as Hexagonal
import ase.io
from ase.lattice.spacegroup import crystal
import copy


class add_strain(object):

    def __init__(self):
        self._delta = 0.02

    def volumetric_strain(self, delta=None):
        # 3c11 + 6C12 #
        if delta is None:
            delta = self._delta
        strain = np.mat([[1 + delta, 0, 0],
                         [0, 1 + delta, 0],
                         [0, 0, 1 + delta]], "float")
        return strain

    def volume_conserving_orthohombic(self, delta=None):
        # 2C11 - 2c12 #
        if delta is None:
            delta = self._delta
        # strain = np.mat([[1 + delta, 0, 0],
        #                  [0, 1 - delta, 0],
        #                  [0, 0, 1 / (1 - delta**2)]], 'float')
        strain = np.mat([[1 / (1 - delta**2), 0., 0.],
                         [0., 1 + delta, 0],
                         [0., 0., 1 - delta]])
        # strain = np.mat([[1 - delta, 0.0, 0.0],
        #                  [0.0, 1. / (1 - delta**2), 0.0],
        #                  [0.0, 0.0, 1 + delta]])
        return strain

    def volume_conserving_monoclinic(self, delta=None):
        # 2c44 #
        if delta is None:
            delta = self._delta

        # strain = np.mat([[1, delta, 0],
        #                  [delta, 1, 0],
        #                  [0, 0, 1 / (1 - delta**2)]], 'float')
        strain = np.mat([[1. / (1 - delta**2), 0.0, 0.0],
                         [0.0, 1.0, delta],
                         [0.0, delta, 1.0]])
        # strain = np.mat([[1.0, 0.0, delta],
        #                  [0.0, 1. / (1 - delta**2), 0.0],
        #                  [delta, 0.0, 1.0]])
        return strain


class bcc(object):

    def set_bcc_primitive_direction(self):
        self._primitive_directions = np.array([[-0.5, 0.5, 0.5],
                                               [0.5, -0.5, 0.5],
                                               [0.5, 0.5, -0.5]])

    def set_bcc_primitive(self, size=(1, 1, 1)):
        atoms = ase.atoms.Atoms(symbols=self.pot["element"],
                                positions=[(0, 0, 0)],
                                info={'unit_cell': 'primitive'},
                                pbc=(1, 1, 1))
        self.set_bcc_primitive_direction()
        cell = self._primitive_directions * self.pot['latbcc']
        atoms.set_cell(cell, scale_atoms=True)
        # atoms = atoms.repeat(size)
        return atoms

    def write_bcc_primitive(self, size=None):
        atoms = self.set_bcc_primitive(size)
        self.write_config_output(atoms)

    def set_bcc_convention(self, directions=None, size=(1, 1, 1)):
        if directions is None:
            directions = self._default_direction
        atoms = Cubic.BodyCenteredCubic(directions=directions,
                                        latticeconstant=self.pot['latbcc'],
                                        size=size,
                                        symbol=self.pot['element'],
                                        pbc=(1, 1, 1))
        atoms.wrap()
        return atoms

    def write_bcc_convention(self, in_direction, size=(1, 1, 1)):
        atoms = self.set_bcc_convention(in_direction, size)
        self.write_config_output(atoms)

    def write_bcc_with_strain(self,
                              delta=None,
                              in_tag='ortho',
                              size=None,
                              write=False):

        if in_tag == 'ortho' or in_tag == 'c12':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono' or in_tag == 'c44':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume' or in_tag == 'c11':
            strain = self.volumetric_strain(delta)

        atoms = self.set_bcc_convention(size)
        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)
        if write is True:
            self.write_config_output(atoms)
        return atoms

    def write_bcc_primitive_with_strain(self,
                                        delta=None,
                                        in_tag=None,
                                        size=None,
                                        strain=None,
                                        write=True):
        if strain is None:
            if in_tag == 'ortho' or in_tag == 'c12':
                strain = self.volume_conserving_orthohombic(delta)
            elif in_tag == 'mono' or in_tag == 'c44':
                strain = self.volume_conserving_monoclinic(delta)
            elif in_tag == 'volume' or in_tag == 'c11':
                strain = self.volumetric_strain(delta)

        atoms = self.set_bcc_primitive(size)
        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)

        if write is True:
            self.write_config_output(atoms)
        return atoms

    def example_of_crystal(self):
        a = 3.3
        Nb = crystal('Nb',
                     [(0, 0, 0)],
                     spacegroup=229,
                     cellpar=[a, a, a, 90, 90, 90],
                     primitive_cell=True)
        print(Nb.get_cell())
        Nb.write('geometry_Bcc.in')


class fcc(object):

    def set_fcc_primitive_direction(self):
        self._primitive_directions = np.array([[0.0, 0.5, 0.5],
                                               [0.5, 0, 0.5],
                                               [0.5, 0.5, 0.0]], 'float')

    def set_fcc_primitive(self, size=(1, 1, 1)):
        atoms = ase.atoms.Atoms(symbols=self.pot["element"],
                                positions=[(0, 0, 0)],
                                info={'unit_cell': 'primitive'},
                                pbc=(1, 1, 1))
        self.set_fcc_primitive_direction()
        cell = self._primitive_directions * self.pot['latfcc']
        atoms.set_cell(cell, scale_atoms=True)
        atoms2 = atoms.repeat(size)
        return atoms2

    def set_fcc_convention(self, directions=None, size=(1, 1, 1)):
        if directions is None:
            l_direction = self._default_direction
        atoms = Cubic.FaceCenteredCubic(directions=l_direction,
                                        latticeconstant=self.pot['latfcc'],
                                        size=size,
                                        symbol=self.pot["element"],
                                        pbc=(1, 1, 1))
        return atoms

    def write_fcc_convention(self, size=None):
        atoms = self.set_fcc_convention(size)
        self.write_config_output(atoms)

    def write_fcc_primitive(self, size=None):
        atoms = self.set_fcc_primitive(size)
        self.write_config_output(atoms)

    def write_fcc_with_strain(self,
                              delta=None,
                              in_tag='ortho',
                              size=None):
        if in_tag == 'ortho':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume':
            strain = self.volumetric_strain(delta)

        atoms = self.set_fcc_convention(size)

        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell
        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)

        self.write_config_output(atoms)

    def write_fcc_primitive_with_strain(self,
                                        delta=None,
                                        in_tag=None,
                                        size=None):
        if in_tag == 'ortho' or in_tag == 'c12':
            strain = self.volume_conserving_orthohombic(delta)
        elif in_tag == 'mono' or in_tag == 'c44':
            strain = self.volume_conserving_monoclinic(delta)
        elif in_tag == 'volume' or in_tag == 'c11':
            strain = self.volumetric_strain(delta)

        atoms = self.set_fcc_primitive(size)

        org_cell = atoms.get_cell()
        new_cell = copy.deepcopy(org_cell)

        org_positions = np.mat(atoms.get_positions())
        org_positions = org_positions * strain

        new_cell = strain * new_cell

        atoms.set_cell(new_cell)
        atoms.set_positions(org_positions)

        self.write_config_output(atoms)


class hcp(object):

    def set_hcp_lattice_constant_ratio(self,
                                       ta=np.sqrt(3) / 2.,
                                       tb=np.sqrt(8. / 3.)):
        self.pot['ahcp'] = self.pot['lattice'] * ta
        self.pot['chcp'] = self.pot['ahcp'] * tb

    def set_hcp_direction(self):
        self.hcp_directions = np.array([[1, 0, 0],
                                        [-0.5, np.sqrt(3) / 2, 0],
                                        [0, 0, np.sqrt(8. / 3)]])

    def set_hcp_convention(self, size=(1, 1, 1)):
        atoms = Hexagonal.HexagonalClosedPacked(
            latticeconstant={'a': self.pot['ahcp'],
                             'c': self.pot['chcp']},
            size=size,
            symbol=self.pot["element"],
            pbc=(1, 1, 1))
        print(atoms.get_cell())
        return atoms

    def write_hcp_convention(self, size=None):
        atoms = self.set_hcp_convention(size=size)
        self.write_config_output(atoms)

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


class gnStructure(bcc, fcc, hcp, add_strain):

    def __init__(self, pot=None):
        add_strain.__init__(self)
        self._default_direction = [[1, 0, 0],
                                   [0, 1, 0],
                                   [0, 0, 1]]

    def mymkdir(self, dirname):
        if not os.path.isdir(dirname):
            os.mkdir(dirname)

    def write_config_with_fix(self, atoms):
        #  a_list = np.array([102, 107, 103, 105, 104]);
        #  b_list = np.array([123, 122, 124, 84,   89]);
        fix_list = np.array([107, 103, 105, 122, 124, 84])
        cell = atoms.get_cell()
        positions = atoms.get_positions()
        print(positions)
        with open("POSCAR", 'w') as fid:
            fid.write("W\n")
            fid.write("1\n")
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[0, 0], cell[0, 1], cell[0, 2]))
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[1, 0], cell[1, 1], cell[1, 2]))
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[2, 0], cell[2, 1], cell[2, 2]))
            fid.write("W\n")
            fid.write("%d\n" % (len(atoms)))
            fid.write("Selective dynamics\n")
            fid.write("Cartesian\n")
            for i in range(len(atoms)):
                if ((i == fix_list).any()):
                    fid.write("%12.9f %12.9f %12.9f   T   T   F\n"
                              % (positions[i, 0],
                                 positions[i, 1],
                                 positions[i, 2]))
                else:
                    fid.write("%12.9f %12.9f %12.9f   T   T   T\n"
                              % (positions[i, 0],
                                 positions[i, 1],
                                 positions[i, 2]))
            fid.close()
        os.system("cp POSCAR POSCAR.vasp")

    def vasp_intro_alloy(self):
        atoms = ase.io.read(filename="POSCAR", format='vasp')
        print(atoms)
        alloy_num = 0
        for atom in atoms:
            if np.random.rand() < 0.1:
                atom.symbol = 'Re'
                alloy_num += 1
        cell = atoms.get_cell()
        with open("POSCAR", 'w') as fid:
            fid.write("W   Re\n")
            fid.write("1\n")
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[0, 0], cell[0, 1], cell[0, 2]))
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[1, 0], cell[1, 1], cell[1, 2]))
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[2, 0], cell[2, 1], cell[2, 2]))
            fid.write("W   Re\n")
            fid.write("%d   %d\n" % ((len(atoms) - alloy_num), alloy_num))
            fid.write("Cartesian\n")
            for atom in atoms:
                if atom.symbol == 'W':
                    fid.write("%12.9f %12.9f %12.9f\n"
                              % (atom.position[0],
                                 atom.position[1],
                                 atom.position[2]))
            for atom in atoms:
                if atom.symbol == 'Re':
                    fid.write("%12.9f %12.9f %12.9f\n"
                              % (atom.position[0],
                                 atom.position[1],
                                 atom.position[2]))
            fid.close()
        os.system("cp POSCAR POSCAR.vasp")

    def cut_plane_for_vacancy(self, atoms):
        cell = atoms.get_cell()
        print(cell[2, 2])
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

        cx = np.dot(in_cell[2, :], vect_x.transpose())

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

    def write_lmp_config_data_charge(self, atoms, filename="lmp_init.txt"):
        pos = atoms.get_positions()
        natm = len(pos)
        cell = atoms.get_cell()
        allsymbos = atoms.get_chemical_symbols()
        symbols = np.unique(allsymbos)
        print(symbols)
        with open(filename, mode="w") as fid:
            fid.write("#lmp data config")
            fid.write("\n")
            fid.write("%d atoms\n" % (natm))
            fid.write("{} atom types\n".format(len(symbols)))
            fid.write("%f\t%f xlo xhi\n" % (0, cell[0, 0]))
            fid.write("%f\t%f ylo yhi\n" % (0, cell[1, 1]))
            fid.write("%f\t%f zlo zhi\n" % (0, cell[2, 2]))
            fid.write("%f  %f  %f xy xz yz\n"
                      % (cell[1, 0], cell[2, 0], cell[2, 1]))
            fid.write("Atoms\n")
            fid.write("\n")
            for i in range(natm):
                for k in range(len(symbols)):
                    if allsymbos[i] == symbols[k]:
                        if k in [0, 1, 2]:
                            fid.write("%d %d %g %12.9f %12.9f %12.9f\n"
                                      % (i + 1, k + 1, 0, pos[i, 0], pos[i, 1], pos[i, 2]))
                        elif k in [3]:
                            fid.write("%d %d %g %12.9f %12.9f %12.9f\n"
                                      % (i + 1, k + 1, 3, pos[i, 0], pos[i, 1], pos[i, 2]))
                        elif k in [4]:
                            fid.write("%d %d %g %12.9f %12.9f %12.9f\n"
                                      % (i + 1, k + 1, -2, pos[i, 0], pos[i, 1], pos[i, 2]))
        fid.close()

    def write_poscar_fix(self, atoms, filename="POSCAR"):
        positions = atoms.get_positions()
        natm = len(positions)
        cell = atoms.get_cell()
        allsymbos = atoms.get_chemical_symbols()
        symbols = np.unique(allsymbos)
        nb = []

        for k in range(len(symbols)):
            cn = 0
            for i in range(natm):
                if (allsymbos[i] == symbols[k]):
                    cn += 1
            nb.append(cn)

        with open(filename, "w") as fid:
            fid.write("poscar\n")
            fid.write("1\n")
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[0, 0], cell[0, 1], cell[0, 2]))
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[1, 0], cell[1, 1], cell[1, 2]))
            fid.write("%12.9f %12.9f %12.9f\n" %
                      (cell[2, 0], cell[2, 1], cell[2, 2]))
            for i in range(len(symbols)):
                fid.write("{} ".format(symbols[i]))
            fid.write("\n")
            for i in range(len(symbols)):
                fid.write("{} ".format(nb[i]))
            fid.write("\n")
            fid.write("Selective dynamics\n")
            fid.write("Cartesian\n")
            for k in range(len(symbols)):
                for i in range(natm):
                    if (allsymbos[i] == symbols[k]):
                        fid.write("%12.9f %12.9f %12.9f F F T\n"
                                  % (positions[i, 0], positions[i, 1], positions[i, 2]))
        fid.close()

    def write_lmp_config_data(self, atoms, filename="lmp_init.txt"):
        positions = atoms.get_positions()
        natm = len(positions)
        cell = atoms.get_cell()
        allsymbos = atoms.get_chemical_symbols()
        symbols = np.unique(allsymbos)
        with open(filename, mode="w") as fid:
            fid.write("#lmp data config")
            fid.write("\n")
            fid.write("%d atoms\n" % (natm))
            fid.write("{} atom types\n".format(len(symbols)))
            fid.write("%f\t%f xlo xhi\n" % (0, cell[0, 0]))
            fid.write("%f\t%f ylo yhi\n" % (0, cell[1, 1]))
            fid.write("%f\t%f zlo zhi\n" % (0, cell[2, 2]))
            fid.write("%f  %f  %f xy xz yz\n"
                      % (cell[1, 0], cell[2, 0], cell[2, 1]))
            fid.write("Atoms\n")
            fid.write("\n")
            for i in range(natm):
                for k in range(len(symbols)):
                    if allsymbos[i] == symbols[k]:
                        fid.write("%d %d %12.9f %12.9f %12.9f\n"
                                  % (i + 1, k + 1, positions[i, 0], positions[i, 1], positions[i, 2]))
        fid.close()

    def write_potfit_config(self, atoms, filename="dummy.config"):
        positions = atoms.get_positions()
        natm = len(positions)
        cell = atoms.get_cell()
        allsymbos = atoms.get_chemical_symbols()
        symbols = np.unique(allsymbos)
        print(symbols)
        with open(filename, mode="w") as fid:
            fid.write("#N {} 1\n".format(natm))
            fid.write("#C Nb\n")
            fid.write("## abc\n")
            fid.write("#X {} {} {}\n".format(
                cell[0, 0], cell[0, 1], cell[0, 2]))
            fid.write("#Y {} {} {}\n".format(
                cell[1, 0], cell[1, 1], cell[1, 2]))
            fid.write("#Z {} {} {}\n".format(
                cell[2, 0], cell[2, 1], cell[2, 2]))
            fid.write("#W 1\n")
            fid.write("#E 1.0\n")
            fid.write("#S 0 0 0 0 0 0\n")
            fid.write("#F\n")
            for i in range(natm):
                for k in range(len(symbols)):
                    if allsymbos[i] == symbols[k]:
                        fid.write("0 %12.9f %12.9f %12.9f 0.0 0.0 0.0\n"
                                  % (positions[i, 0], positions[i, 1], positions[i, 2]))
        fid.close()

    def write_lmp_coords(self, atoms, filename="lmp_coord"):
        positions = atoms.get_positions()
        natm = len(positions)
        with open(filename, mode="w") as fid:
            fid.write("%d\n" % (natm))
            for i in range(natm):
                fid.write("%d  %f  %f  %f \n"
                          % (i + 1, positions[i, 0], positions[i, 1], positions[i, 2]))
            fid.close()
