#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-13 23:24:50


from numpy import sqrt
import numpy as np
import pandas as pd

Mg = {'pois': 0.29,
      'shearModulus': 17,
      'cOVERa': 1.6236}

c = Mg['cOVERa']

plnsb1 = {
    'Ba': np.array([0, 0, 1]),
    'Pr': np.array([sqrt(3), 1, 0]),
    'PyI': np.array([sqrt(3) * c, c, -sqrt(3)]),
    'PyII': np.array([sqrt(3), c, sqrt(3)])
}

plnsb2 = {
    'Ba': np.array([0, 0, 1]),
    'Pr': np.array([0, 1, 0]),
    'PyI': np.array([0, 2 * c, sqrt(3)]),
    'PyII': np.array([0, -2 * c, sqrt(3)])
}

plnsb3 = {
    'Ba': np.array([0, 0, 1]),
    'Pr': np.array([sqrt(3), -1, 0]),
    'PyI': np.array([sqrt(3) * c, -c, sqrt(3)]),
    'PyII': np.array([-sqrt(3) * c, c, sqrt(3)])
}

plnsb4 = {
    'Pr': np.array([sqrt(3), 1., 0]),
    'PyI': np.array([sqrt(3) * c, -c, sqrt(3)]),
    'PyII': np.array([0, -2 * c, sqrt(3)]),
    'sP': np.array([c, -sqrt(3) * c, 2])}


hcpburgs = {'b1': [np.array([-1. / 2., sqrt(3) / 2., 0.]),
                   plnsb1],
            'b2': [np.array([1., 0., 0.]),
                   plnsb2],
            'b3': [np.array([1. / 2., sqrt(3.) / 2., 0.0]),
                   plnsb3],
            'b4': [np.array([-1. / 2., sqrt(3.) / 2., c]),
                   plnsb4]}

# df = pd.DataFrame(burgs1, index=('burger'), columns=burgcolm)
burgdf = pd.DataFrame(hcpburgs)


class dd_dat:
    domid = 0
    cell = None
    nnodes = None
    datfilev = 5
    precn = None


class arm(object):
    armnode = None
    burg = None
    plane = None


class node(object):
    domid = None
    nodeid = None
    pos = None
    narm = None
    const = 0
    arml = None
    armr = None


class prec(object):
    precid = None
    coords = None
    dimaxi = None
    rotate = None
    strain = None


precfile_header = """# Data for an inclusion consists of the following
# items separated # white space
#
#     Inclusion ID:  Unique integer identifying the inclusion regardless
#                    of its position in the simulation or the domain
#                    encompassing it.  (Should be sequential values)
#     Position:      3 values specifying the XYZ coordinates of
#                    the center of the inclusion
#     Semi-principal axes: 3 values in units of b
#     Euler angles:  3 values specifying the three Euler angles. 1st and 3rd in [-180, 180]. 2nd in [0, 180].
#     Strain field:  Strain field for the particle.  Strain field is a
#                    symetric matrix so only six components are specified
#                    in the order: S[0][0], S[1][1], S[2][2], S[1][2],
#                    S[0][2], S[0][1]
"""
