# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-09 21:23:40
# @Last Modified by:   chaomy
# @Last Modified time: 2017-12-10 14:06:15

import numpy as np

class gn_dd_cell(object):

    def set_cell(self, side=3e3):
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
