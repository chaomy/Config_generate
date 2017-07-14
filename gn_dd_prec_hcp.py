#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-10 08:37:35
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-13 20:27:31

import numpy as np
from numpy import cos, sin, pi
import gn_dd_data_dat as dddat


class gn_dd_prec(object):

    def __init__(self):
        self.cell = np.ndarray([3, 2])
        self.precs = None
        return

    def inplane_hcp_beta1_prec(self):
        self.ddata.precn = 30
        self.set_cell()
        self.set_prec()
        fid = self.write_precip_header()
        self.write_precip_data(fid)
        return

    def set_prec_strain(self):
        strain = np.zeros(6)
        return strain

    def set_prec_rotate(self, angle=(0., 0., 0.)):
        # check  http://mathworld.wolfram.com/EulerAngles.html
        (phi, theta, psi) = angle[0], angle[1], angle[2]
        D = np.mat([[cos(phi), sin(phi), 0.],
                    [-sin(phi), cos(phi), 0.],
                    [0., 0., 1.]])
        C = np.mat([[1., 0., 0.],
                    [0., cos(theta), sin(theta)],
                    [0., -sin(theta), cos(theta)]])
        B = np.mat([[cos(psi), sin(psi), 0.],
                    [-sin(psi), cos(psi), 0.],
                    [0., 0., 1.]])
        A = B * C * D
        return A

    def set_prec_size(self):
        lens = np.array([50, 150, 150]) * np.random.rand(3)
        size = np.array([200, 800, 800]) + lens
        return size

    def set_prec_coords(self):
        cell = self.ddata.cell
        pos = np.zeros(3)
        size = cell[:, 1] - cell[:, 0]
        pos = cell[:, 0] + size * np.random.rand(3)
        # along the (x, y, 0) plane
        pos[2] = 0.0
        return pos

    def set_prec(self):
        self.precs = []
        for i in range(self.ddata.precn):
            prec = dddat.prec()
            prec.precid = i + 1
            prec.coords = self.set_prec_coords()
            prec.dimaxi = self.set_prec_size()
            if i % 3 == 0:
                prec.rotate = self.set_prec_rotate((pi / 3., 0., 0.))
            elif i % 3 == 1:
                prec.rotate = self.set_prec_rotate((-pi / 3., 0., 0.))
            elif i % 3 == 2:
                prec.rotate = self.set_prec_rotate((0., 0., 0.))
            prec.strain = self.set_prec_strain()
            print "coord", prec.coords
            print "size", prec.dimaxi
            self.precs.append(prec)
        return

    def write_precip_header(self, fname='prect.dat'):
        fid = open(fname, 'w')
        fid.write(dddat.precfile_header)
        return fid

    def write_precip_data(self, fid):
        strformat = '{} ' + '{:4.2f} ' * (3 * 6) + '\n'
        prec = self.precs[1]
        for prec in self.precs:
            line = strformat.format(
                prec.precid,
                prec.coords[0], prec.coords[1], prec.coords[2],
                prec.dimaxi[0], prec.dimaxi[1], prec.dimaxi[2],
                prec.rotate[0, 0], prec.rotate[0, 1],
                prec.rotate[0, 2], prec.rotate[1, 0],
                prec.rotate[1, 1], prec.rotate[1, 2],
                prec.strain[0], prec.strain[1],
                prec.strain[2], prec.strain[3],
                prec.strain[4], prec.strain[5])
            fid.write(line)
        fid.close()
        return
