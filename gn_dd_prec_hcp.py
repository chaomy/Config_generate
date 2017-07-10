#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-10 08:37:35
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-10 09:40:09

import numpy as np
import gn_dd_data_dat as dddat


class gn_dd_prec():

    def __init__(self):
        self.cell = np.ndarray([3, 2])
        self.precs = None
        return

    def inplane_hcp_beta1_prec(self):
        self.set_cell()
        self.ddata.precn = 10
        self.set_prec_pos()
        return

    def set_prec_pos(self):
        self.precs = []
        for i in range(self.ddata.precn):
            prec = ddata.prec
            self.precs.append(dddat.prec)
        print self.precs
        return

    def write_precip_header(self):

        return

    def write_precip_data(self):

        return
