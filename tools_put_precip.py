#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-25 14:28:58
# @Last Modified by:   yangchaoming
# @Last Modified time: 2017-07-10 08:09:43


import math
import numpy as np


class paradis_precip_tools(object):

    def __init__(self):
        self.r = 300
        self.p1 = np.array([0., 0., 0.])
        self.p2 = np.array([math.sqrt(3.) / 2. * self.r,
                            0.5 * self.r,
                            0.0])

        self.b = np.array([-0.5, math.sqrt(3) / 2., 0.0])
        self.t = self.p2 - self.p1
        self.t /= np.linalg.norm(self.t)

        self.n = np.array([0.0, 0.0, 1.0])
        return

    def tools_find_precip_position(self):
        r = 300

        pA = np.array([0., 0., 0.])
        pB = np.array([math.sqrt(3.) / 2. * r, 0.5 * r, 0.0])
        pn = np.array([0., 0., 1.0])

        t = pB - pA

        norm_v = np.cross(t, pn)
        norm_v /= np.linalg.norm(norm_v)
        norm_v *= 500

        cent = 0.5 * (pB + pA)
        pos = cent + norm_v
        print(pos)
        return

    def tools_rotate_shear(self):
        shear_xy = np.mat([[0.0, -0.3, 0.0],
                           [-0.3, 0.0, 0.0],
                           [0.0, 0.0, 0.0]])

        b = self.b
        t = self.t
        #  n = self.n;
        n = np.cross(b, t)

        new_cord = np.mat([[b[0], b[1], b[2]],
                           [t[0], t[1], t[2]],
                           [n[0], n[1], n[2]]
                           ])

        print(new_cord)
        new_shear = new_cord.transpose() * -shear_xy * new_cord
        print(1e9 * new_shear)
        return

if __name__ == '__main__':
    job = paradis_precip_tools()
    #  job.tools_rotate_shear();
    job.tools_find_precip_position()
