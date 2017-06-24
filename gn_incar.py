#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./cal_vasp_cij.py
#
###################################################################
#
# Purpose :  generate the incar for vasp
# write_incar
#
# Creation Date :
# Last Modified : Sat Apr  1 23:15:50 2017
# Created By    : Chaoming Yang
#
###################################################################


class gn_incar(object):

    def __init__(self, in_type=None):
        print in_type
        if in_type is not None:
            self.in_type = in_type
        else:
            self.in_type = "dft"

        self.nsw = 1000
        self._ab_md_temp = 1000
        self._accuracy = 1e-5
        self._encut = 550
        self._enaug = 550
        self._npar = 4
        self._ibrion = 2
        return

    def set_nsw(self, nsw):
        self.nsw = nsw
        return

    def set_ibrion(self, ibrion):
        self._ibrion = ibrion
        return

    def set_accuracy(self, accuracy):
        self._accuracy = accuracy
        return

    def set_ab_md_temp(self, temp):
        self._ab_md_temp = temp
        return

    def set_encut(self, in_ecut, in_eaug):
        self._encut = in_ecut
        self._enaug = in_eaug
        return

    def set_incar_type(self, in_type):
        self.in_type = in_type
        return

    def set_npar(self, npar):
        self._npar = npar
        return

    def write_incar(self, fname="INCAR"):
        if self.in_type == 'dftunrelax' or self.in_type == 'scf':
            with open(fname, 'w') as fid:
                fid.write("""
PREC   = A
IBRION =-1
ISMEAR =-5
ISIF   = 2

ENCUT  = 500.00000   # 1.5*(Max ENMAX)
EDIFF  = %1.0e       ## 1e-6
NPAR   = %d

LCHARG = .FALSE.
LWAVE  = .FALSE.
LELF   = .FALSE.
LVTOT  = .FALSE.
""" % (self._accuracy, self._npar))
                fid.close()

        elif self.in_type == 'dftrelax' or self.in_type == 'relax':
            with open(fname, 'w') as fid:
                fid.write("""
NSW    = %d
IBRION = 2
ISIF   = 2
ISMEAR = 1

LELF   = .FALSE.
LCHARG = .FALSE.
LWAVE  = .FALSE.
LVTOT  = .FALSE.
NPAR   = %d
""" % (self.nsw, self._accuracy, self._npar))
                fid.close()

        elif self.in_type == 'ab_md': 
            with open(fname, 'w') as fid:
                fid.write("""
PREC   = Med
NSW    = {} 
IBRION = 0      
POTIM  = 0.5    
ISIF   = 2     
TEBEG  = {}
TEEND  = {}
SMASS  = 0

ISMEAR = 1
EDIFF  = 1e-2
NPAR   = 4

LELF   = .FALSE.
LCHARG = .FALSE.
LWAVE  = .FALSE.
LVTOT  = .FALSE.
                         """ % (self.nsw,
                                self._ab_md_temp,
                                self._ab_md_temp))
                fid.close()

        elif self.in_type == 'phonon':
            with open(fname, 'w') as fid:
                fid.write("""
NWRITE = 2
PREC   = A

NSW    = 1
IBRION = 8
ISIF   = 2
TEBEG  = 300.

ISMEAR = 1
LELF   = T
LVTOT  = T

IALGO  = 38
LREAL  = .FALSE.
ISPIN  = 1
ENCUT  = 400
EDIFF  = 1e-8
ADDGRID = .TRUE.
LWAVE  = .FALSE.
LCHARG = .FALSE.
LVTOT  = .FALSE.
NELM   = 400
NBANDS = 118
""")
                fid.close()
        return


if __name__ == '__main__':
    from optparse import OptionParser
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-t", "--type",
                      action="store",
                      type="string",
                      dest="type")
    (options, args) = parser.parse_args()

    M = gn_incar(options.type)
    M.write_incar()
