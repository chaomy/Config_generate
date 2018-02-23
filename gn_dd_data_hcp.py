#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-02-19 02:06:01


import numpy as np
import gn_dd_data_dat as dddat
import gn_dd_load
import gn_dd_ctrl
import gn_dd_prec
import gn_dd_cell


class gn_dd_data_hcp(gn_dd_prec.gn_dd_prec,
                     gn_dd_load.gn_dd_load,
                     gn_dd_ctrl.gn_dd_ctrl,
                     gn_dd_cell.gn_dd_cell):

    def __init__(self):
        self.job = 'paradis'
        gn_dd_prec.gn_dd_prec.__init__(self)
        gn_dd_load.gn_dd_load.__init__(self)
        gn_dd_ctrl.gn_dd_ctrl.__init__(self)
        gn_dd_cell.gn_dd_cell.__init__(self)
        self.ddata = dddat.dd_dat
        self.ddata.nnodes = 6
        self.burgs = dddat.hcpslip
        self.bkey = 'b1'
        self.pln = 'Ba'
        self.datafile = '{}.data'.format(self.job)

    def gn_hcp_straight_dis(self):
        nnodes = self.ddata.nnodes
        bkey = self.bkey
        self.set_cell()

        # write
        length = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0])
        disl = 0.3 * length
        delta = disl / (nnodes - 1)

        idr = range(nnodes)
        idl = range(nnodes)
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())

        nlist = []
        plane = self.burgs[bkey]['norms']['Ba']
        for i, lid, rid in zip(range(nnodes), idl, idr):
            node = dddat.node()
            node.domid = self.ddata.domid
            node.nodeid = i
            node.pos = np.zeros(3)

            node.pos[0] = \
                -0.5 * disl + i * delta
            node.narm = 2
            node.arml = self.set_arm(lid,
                                     self.burgs[bkey]['b'], plane)
            node.armr = self.set_arm(rid,
                                     -self.burgs[bkey]['b'], plane)
            nlist.append(node)
            if i in [0, nnodes - 1]:  # add constrain
                node.const = 7
                node.narm = 1
        return nlist

    def write_hcp_straight_data(self):
        nlist = self.gn_hcp_straight_dis()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data(nlist, fid)

    def write_hcp_orawan_data(self):
        self.set_fname('hcp')
        self.write_ctrl_file()  # write ctrl file

        # write data file
        nlist = self.gn_hcp_straight_dis()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data_end(nlist, fid)

        self.inplane_hcp_beta1_prec()  # write preicpitates

    def write_hcp_tensile_data(self):
        self.set_fname('mgp016r1e5')
        # write ctrl file
        self.write_ctrl_file(ltype='strain')
        self.hcp_beta1_prec()
