#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-31 16:35:33


import numpy as np
import gn_dd_data_dat as dddat
import gn_dd_load
import gn_dd_ctrl
import gn_dd_prec


class gn_dd_data_hcp(gn_dd_prec.gn_dd_prec,
                     gn_dd_load.gn_dd_load,
                     gn_dd_ctrl.gn_dd_ctrl):

    def __init__(self):
        self.job = 'paradis'
        gn_dd_prec.gn_dd_prec.__init__(self)
        gn_dd_load.gn_dd_load.__init__(self)
        gn_dd_ctrl.gn_dd_ctrl.__init__(self)
        self.ddata = dddat.dd_dat
        self.ddata.nnodes = 4
        self.burgs = dddat.hcpslip
        self.bkey = 'b1'
        self.pln = 'Ba'
        self.datafile = '{}.data'.format(self.job)
        return

    def set_domid(self, domid=0):
        self.ddata.domid = domid
        return

    def set_cell(self, side=9e3):
        cell = np.ndarray([3, 2])
        cell[:, 0] = np.ones(3) * -side
        cell[:, 1] = np.ones(3) * side
        self.ddata.cell = cell
        self.cal_cell_vol()
        return

    def cal_cell_vol(self):
        cell = self.ddata.cell
        self.ddata.cellvol = np.prod(cell[:, 1] - cell[:, 0])
        return

    def set_arm(self, armid, burg, plane):
        arm = dddat.arm()
        arm.armnode = armid
        arm.burg = burg
        arm.plane = plane
        return arm

    def gn_hcp_straight_dis(self):
        nnodes = self.ddata.nnodes
        bkey = self.bkey
        self.set_cell()
        length = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0])
        disl = 0.2 * length
        delta = disl / (nnodes - 1)
        idr = range(nnodes)
        idl = range(nnodes)
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())
        nlist = []
        plane = self.burgs[bkey]['norms']['PyI']
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
            if i in [0, nnodes - 1]:
                node.const = 7
                node.narm = 1
        return nlist

    def write_data_head(self):
        cell = self.ddata.cell
        nnodes = self.ddata.nnodes
        fid = open(self.datafile, 'w')
        fid.write(""" dataFileVersion =   5
numFileSegments =   1
minCoordinates = [
  {}
  {}
  {}
  ]
maxCoordinates = [
  {}
  {}
  {}
  ]
nodeCount = {}
dataDecompType =   2
dataDecompGeometry = [
  1
  1
  1
  ]
 """.format(cell[0, 0], cell[1, 0],
            cell[2, 0], cell[0, 1],
            cell[1, 1], cell[2, 1], nnodes))
        return fid

    def write_domain_data(self, fid):
        cell = self.ddata.cell
        fid.write("""domainDecomposition =
# Dom_ID  Minimum XYZ bounds   Maximum XYZ bounds
  {}  {}  {}  {}  {}  {}  {}
 """.format(self.ddata.domid,
            cell[0, 0], cell[1, 0], cell[2, 0],
            cell[0, 1], cell[1, 1], cell[2, 1]))
        return fid

    def write_arm_data(self, arm, fid):
        fid.write("  {},{} {:8.7f} {:8.7f} {:8.7f}\n".format(self.ddata.domid,
                                                             arm.armnode,
                                                             arm.burg[0],
                                                             arm.burg[1],
                                                             arm.burg[2]))
        fid.write("  	  {:8.7f} {:8.7f} {:8.7f}\n".format(arm.plane[0],
                                                          arm.plane[1],
                                                          arm.plane[2]))
        return fid

    def write_nodal_data(self, nlist, fid):
        fid.write("nodalData = \n")
        fid.write("# Primary lines: node_tag, x, y, z, num_arms, constraint \n")
        fid.write("# Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n")
        for node in nlist:
            fid.write("{},{} {} {} {} {} {}\n".format(node.domid,
                                                      node.nodeid, node.pos[0],
                                                      node.pos[1], node.pos[2],
                                                      node.narm, node.const))
            fid = self.write_arm_data(node.arml, fid)
            fid = self.write_arm_data(node.armr, fid)
        return fid

    def write_nodal_data_end(self, nlist, fid):
        fid.write("nodalData = \n")
        fid.write("# Primary lines: node_tag, x, y, z, num_arms, constraint \n")
        fid.write("# Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n")
        for node in nlist:
            fid.write("{},{} {} {} {} {} {}\n".format(node.domid,
                                                      node.nodeid, node.pos[0],
                                                      node.pos[1], node.pos[2],
                                                      node.narm, node.const))
            if node.nodeid in [0]:
                fid = self.write_arm_data(node.armr, fid)
            elif node.nodeid in [self.ddata.nnodes - 1]:
                fid = self.write_arm_data(node.arml, fid)
            else:
                fid = self.write_arm_data(node.arml, fid)
                fid = self.write_arm_data(node.armr, fid)
        return fid

    def write_hcp_straight_data(self):
        nlist = self.gn_hcp_straight_dis()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data(nlist, fid)
        return

    def write_hcp_orawan_data(self):
        self.set_fname('hcp')

        # write ctrl file
        self.write_ctrl_file()

        # write data file
        nlist = self.gn_hcp_straight_dis()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data_end(nlist, fid)

        # write preicpitates
        self.inplane_hcp_beta1_prec()
        return

    def write_hcp_tensile_data(self):
        self.set_fname('mgprec016')
        # write ctrl file
        self.write_ctrl_file(ltype='strain')
        self.hcp_beta1_prec()
        return
