#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-14 15:05:33


import numpy as np
from . import gn_dd_data_dat as dddat
from . import gn_dd_load
from . import gn_dd_ctrl
from . import gn_dd_prec
from . import gn_dd_cell


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
        self.ndis = 7
        self.ddata.nnodes = 20
        self.ddata.totalnodes = self.ddata.nnodes * self.ndis
        self.burgs = dddat.hcpslip
        self.bkey = 'b1'
        self.pln = 'Ba'
        self.datafile = '{}.data'.format(self.job)

    def gn_hcp_straight_dis_periodic(self):
        nnodes = self.ddata.nnodes
        bkey = self.bkey
        self.set_cell()

        length = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0])
        disl = length
        delta = disl / (nnodes)

        idr = list(range(nnodes))
        idl = list(range(nnodes))
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())

        nlist = []
        plane = self.burgs[bkey]['norms']['Ba']
        for i, lid, rid in zip(list(range(nnodes)), idl, idr):
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
        return nlist

    def gn_hcp_straight_dis_many(self):
        nnodes = self.ddata.nnodes
        bkey = self.bkey
        self.set_cell()

        length = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0])
        disl = length
        delta = disl / (nnodes - 1)

        idr = list(range(nnodes))
        idl = list(range(nnodes))
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())

        nlist = []
        plane = self.burgs[bkey]['norms']['Ba']

        for k in range(self.ndis):
            for i, lid, rid in zip(list(range(nnodes)), idl, idr):
                node = dddat.node()
                node.domid = self.ddata.domid
                node.nodeid = i + nnodes * k
                node.pos = np.zeros(3)
                node.pos[0] = \
                    -0.5 * disl + i * delta
                node.pos[1] = 50 * k + 2800
                node.narm = 2
                node.arml = self.set_arm(lid + nnodes * k,
                                         self.burgs[bkey]['b'], plane)
                node.armr = self.set_arm(rid + nnodes * k,
                                         -self.burgs[bkey]['b'], plane)
                nlist.append(node)
        return nlist

    def gn_hcp_straight_dis_fixends(self):
        nnodes = self.ddata.nnodes
        bkey = self.bkey
        self.set_cell()

        length = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0])
        disl = 0.4 * length
        delta = disl / (nnodes - 1)

        idr = list(range(nnodes))
        idl = list(range(nnodes))
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())

        nlist = []
        plane = self.burgs[bkey]['norms']['Ba']
        for i, lid, rid in zip(list(range(nnodes)), idl, idr):
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
        self.set_fname('hcp')
        self.write_ctrl_file()
        # nlist = self.gn_hcp_straight_dis_periodic()
        nlist = self.gn_hcp_straight_dis_many()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data(nlist, fid)

    def write_hcp_orawan_data(self):
        self.set_fname('hcp')
        self.write_ctrl_file()  # write ctrl file
        nlist = self.gn_hcp_straight_dis_fixends()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data_end(nlist, fid)
        # self.inplane_hcp_beta1_prec()  # write preicpitates

    def write_hcp_tensile_data(self):
        self.set_fname('mgpr1e5')
        self.write_ctrl_file(ltype='strain')
        self.hcp_beta1_prec()

# paradisrepart -infile paradis.data -cells 4,4,4 -domains 4,4,4 -outfile
# output.data -domainType  2
