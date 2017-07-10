#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-10 00:46:03


import numpy as np
from itertools import cycle
import gn_dd_data_dat as dddat


class gn_dd_data_hcp(object):

    def __init__(self):
        self.ddata = dddat.dd_dat
        self.ddata.nnodes = 6
        return

    def set_domid(self, domid=0):
        self.ddata.domid = domid
        return

    def set_cell(self):
        cell = np.ndarray([3, 2])
        cell[:, 0] = np.ones(3) * -1e3
        cell[:, 1] = np.ones(3) * 1e3
        self.ddata.cell = cell
        return

    def set_arm(self, armid, burg, plane):
        arm = dddat.arm()
        arm.armnode = armid
        arm.burg = burg
        arm.plane = plane
        return arm

    def gn_hcp_straight_dis(self):
        nnodes = self.ddata.nnodes
        drv.set_cell()
        delta = (self.ddata.cell[0, 1] - self.ddata.cell[0, 0]) / nnodes
        idr = range(nnodes)
        idl = range(nnodes)
        idr.append(idr.pop(0))
        idl.insert(0, idl.pop())
        burgl = np.array([1. / 2., np.sqrt(3.) / 2., 0.0])
        burgr = -np.array([1. / 2., np.sqrt(3.) / 2., 0.0])
        plane = np.array([0, 0, 1])
        nlist = []

        for i, lid, rid in zip(range(nnodes), idl, idr):
            node = dddat.node()
            node.domid = self.ddata.domid
            node.nodeid = i
            node.pos = np.zeros(3)
            node.pos[0] = self.ddata.cell[0, 0] + i * delta
            node.narm = 2
            node.arml = self.set_arm(lid, burgl, plane)
            node.armr = self.set_arm(rid, burgr, plane)
            nlist.append(node)
        return nlist

    def write_data_head(self):
        cell = self.ddata.cell
        nnodes = self.ddata.nnodes
        fid = open('paradis.data', 'w')
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

    def write_hcp_straight_data(self):
        nlist = self.gn_hcp_straight_dis()
        fid = self.write_data_head()
        fid = self.write_domain_data(fid)
        fid = self.write_nodal_data(nlist, fid)
        return


if __name__ == '__main__':
    drv = gn_dd_data_hcp()
    drv.write_hcp_straight_data()