#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-10 08:37:35
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-28 16:42:58

import numpy as np
from numpy import cos, sin, pi
import gn_dd_data_dat as dddat


class gn_dd_prec(object):

    def __init__(self):
        self.precs = None
        self.precfile = '{}.dat'.format(self.job)
        return

    def set_num_prec(self, volfrac=0.16):
        vol0 = volfrac * self.ddata.cellvol
        print vol0
        volsum = 0.0
        cnt = 0
        while True:
            volsum += np.prod(self.set_prec_size())
            print volsum
            cnt += 1
            if (volsum > vol0):
                break
        print "num precs", cnt
        return

    # precipiates in signle basal plane
    def inplane_hcp_beta1_prec(self):
        self.ddata.precn = 66
        self.set_cell()
        self.set_inplane_prec()
        fid = self.write_precip_header()
        self.write_precip_data(fid)
        return

    # precipiates in 3d distribution
    def hcp_beta1_prec(self):
        self.set_cell()
        self.set_num_prec()
        self.set_3d_prec()
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
        print "rotata", phi, A
        return A

    def set_prec_size(self):
        lens = np.array([40, 2, 40]) * np.random.rand(3)
        size = np.array([603, 30, 603]) + lens
        return size

    def set_prec_3d_coords(self):
        cell = self.ddata.cell
        pos = np.zeros(3)
        size = cell[:, 1] - cell[:, 0]
        pos = cell[:, 0] + size * np.random.rand(3)
        return pos

    def set_prec_inplane_coords(self):
        cell = self.ddata.cell
        pos = np.zeros(3)
        size = cell[:, 1] - cell[:, 0]
        pos = cell[:, 0] + size * np.random.rand(3)
        # along the (x, y, 0) plane
        pos[2] = 0.0
        return pos

    def set_3d_prec(self, volfrac=0.016):
        self.precs = []
        vol0 = volfrac * self.ddata.cellvol
        volsum = 0.0
        cnt = 0
        while True:
            cnt += 1
            prec = self.set_beta1(cnt)
            volsum += np.prod(prec.dimaxi)
            self.precs.append(prec)
            if all([(volsum > vol0), (cnt % 3 == 0)]) is True:
                break
        self.ddata.precn = cnt
        return

    def set_beta1(self, cnt):
        prec = dddat.prec()
        prec.precid = cnt
        prec.coords = self.set_prec_3d_coords()
        prec.dimaxi = self.set_prec_size()
        if cnt % 3 == 0:
            prec.rotate = self.set_prec_rotate((pi / 3., 0., 0.))
        elif cnt % 3 == 1:
            prec.rotate = self.set_prec_rotate((pi * 2. / 3., 0., 0.))
        elif cnt % 3 == 2:
            prec.rotate = self.set_prec_rotate((0., 0., 0.))
        prec.strain = self.set_prec_strain()
        return prec

    def set_inplane_prec(self):
        self.precs = []
        for i in range(self.ddata.precn):
            prec = dddat.prec()
            prec.precid = i + 1
            prec.coords = self.set_prec_inplane_coords()
            prec.dimaxi = self.set_prec_size()
            if i % 3 == 0:
                prec.rotate = self.set_prec_rotate((pi / 3., 0., 0.))
            elif i % 3 == 1:
                prec.rotate = self.set_prec_rotate((pi * 2. / 3., 0., 0.))
            elif i % 3 == 2:
                prec.rotate = self.set_prec_rotate((0., 0., 0.))
            prec.strain = self.set_prec_strain()
            print "coord", prec.coords
            print "size", prec.dimaxi
            self.precs.append(prec)
        return

    def write_precip_header(self):
        fid = open(self.precfile, 'w')
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
