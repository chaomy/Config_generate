# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-09 21:16:29
# @Last Modified by:   chaomy
# @Last Modified time: 2018-12-14 21:18:53

import numpy as np
import gn_dd_data_dat as dddat
import gn_dd_load
import gn_dd_ctrl
import gn_dd_prec
import gn_dd_cell


class gn_dd_data_bcc(gn_dd_prec.gn_dd_prec,
                     gn_dd_load.gn_dd_load,
                     gn_dd_ctrl.gn_dd_ctrl,
                     gn_dd_cell.gn_dd_cell):

    def __init__(self):
        self.job = 'paradis'
        gn_dd_cell.gn_dd_cell.__init__(self)
        gn_dd_prec.gn_dd_prec.__init__(self)
        gn_dd_load.gn_dd_load.__init__(self)
        gn_dd_ctrl.gn_dd_ctrl.__init__(self)
        self.ddata = dddat.dd_dat
        self.ddata.nnodes = 5
        self.burgs = dddat.bccslip
        self.bkey = 'b1'
        self.pln = 'n1'
        self.datafile = '{}.data'.format(self.job)

    def gn_bcc_straight_dis(self):
        nnodes = self.ddata.nnodes
        bkey = self.bkey
        self.set_cell()

        idr = list(range(nnodes))
        idl = list(range(nnodes))
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())

        lth = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0])
        dlt = lth / (nnodes)  # due to periodicity

        # make a [1, 1, 1] dislocation line
        inip = np.array([self.ddata.cell[0, 0],
                         self.ddata.cell[1, 0],
                         self.ddata.cell[2, 0]])
        dltv = np.array([dlt, dlt, dlt])

        print(idl)
        print(idr)
        nlist = []
        plane = self.burgs[bkey]['norms']['n1']

        for i, lid, rid in zip(list(range(nnodes)), idl, idr):
            node = dddat.node()
            node.domid = self.ddata.domid
            node.nodeid = i
            node.pos = inip + i * dltv
            node.narm = 2
            node.arml = self.set_arm(lid, self.burgs[bkey]['b'], plane)
            node.armr = self.set_arm(rid, -self.burgs[bkey]['b'], plane)
            nlist.append(node)
        return nlist

    def write_straight_screw_data(self):
        nlist = self.gn_bcc_straight_dis()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data(nlist, fid)
        self.write_ctrl_file(ltype='stress')
        # self.gn_bcc_single_sphere()
        self.gn_fcc_inplane_sphere()
