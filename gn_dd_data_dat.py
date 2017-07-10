#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-10 00:28:21


class dd_dat:
    domid = 0
    cell = None
    nnodes = None
    datfilev = 5


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
