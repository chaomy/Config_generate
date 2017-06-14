#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ./gn_kpoints.py
#
###################################################################
#
# Purpose :   generate Kpoints for vasp
#
# Creation Date :
# Last Modified : Sat Apr  1 23:15:41 2017
# Created By    : Chaoming Yang
#
###################################################################


class gn_kpoints(object):
    def __init__(self,
                 in_type='gamma',
                 structure='bcc',
                 kpoints=[31, 31, 31]):

        self.in_type = in_type
        self.structure = structure
        kpoints = [31, 31, 31]
        self.kpoints = kpoints
        return

    def set_kpoints(self, k_point):
        self.kpoints = [k_point, k_point, k_point]
        return

    def set_diff_kpoints(self, kpoint_list):
        self.kpoints = kpoint_list
        return

    def set_intype(self, intype):
        self.in_type = intype
        return

    def write_kpoints(self):
        if self.in_type == 'band':
            if self.structure == "bcc":
                with open("KPOINTS", 'w') as fid:
                    fid.write("""#%s
11
Line-mode
Cartesian
0.0 0.0 0.0 ! G
0.0 1.0 0.0 ! H

0.0 1.0 0.0 ! H
0.5 0.5 0.0 ! N

0.5 0.5 0.0 ! N
0.0 0.0 0.0 ! G

0.0 0.0 0.0 ! G
0.5 0.5 0.5 ! P
""" % (self.structure))
                    fid.close()
            if self.structure == "fcc":
                with open("KPOINTS", 'w') as fid:
                    fid.write("""#%s
11
Line-mode
Cartesian
0.5 0.5 0.5 ! L
0.0 0.0 0.0 ! G

0.0 0.0 0.0 ! G
0.0 0.0 1.0 ! X

0.0 0.0 1.0 ! X
0.5 0.0 1.0 ! W

0.5 0.0 1.0 ! W
0.75 0.0 0.75 ! K

0.75 0.0 0.75 ! K
0.0 0.0 0.0 ! G
""" % (self.structure))
                    fid.close()
            if self.structure == "hcp":
                with open("KPOINTS", 'w') as fid:
                    fid.write("""#%s
11
Line-mode
Reciprocal
0.0 0.0 0.0 ! G
0.0 0.0 0.5 ! A

0.0 0.0 0.5 ! A
0.33333 0.33333 0.5 ! H

0.33333 0.33333 0.5 ! H
0.33333 0.33333 0.0 ! K

0.33333 0.33333 0.0 ! K
0.0 0.0 0.0 ! G

0.0 0.0 0.0 ! G
0.0 0.5 0.0 ! M

0.0 0.5 0.0 ! M
0.0 0.5 0.5 ! L

0.0 0.5 0.5 ! L
0.0 0.0 0.0 ! G
 """ % (self.structure))

        if self.in_type != 'band':
            if self.in_type.lower() == 'gamma':
                if self.structure == 'bcc':
                    with open("KPOINTS", 'w') as fid:
                        fid.write("""Automatic mesh
0
Gamma
%d  %d  %d
0  0  0
                            """ % (self.kpoints[0], self.kpoints[1],
                                   self.kpoints[2]))
        else:
            with open("KPOINTS", 'w') as fid:
                fid.write("""Automatic mesh
0
Monkhorst-Pack
%d  %d  %d
0  0  0
""" % (self.kpoints[0], self.kpoints[1], self.kpoints[2]))


if __name__ == '__main__':
    from optparse import OptionParser
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-k",
                      "--KPOINTS",
                      dest="KPOINTS",
                      type="int",
                      help="KPOINTS number",
                      default="33")
    parser.add_option("-t",
                      "--type",
                      action="store",
                      type="string",
                      dest="type",
                      default="lattice")
    parser.add_option(
        "-s", "--structure", action="store", dest='structure', default="bcc")
    (options, args) = parser.parse_args()
    M = gn_kpoints(options.type, options.structure, options.KPOINTS)
    M.write_kpoints()
