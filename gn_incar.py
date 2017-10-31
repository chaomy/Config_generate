# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-09 16:05:24
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-20 14:46:51


class gn_incar(object):

    def __init__(self, in_type='scf'):
        self.in_type = in_type
        self.nsw = 1000
        self._ab_md_temp = 1000
        self._accuracy = 1e-5
        self._encut = 550
        self._enaug = 550
        self._npar = 4
        self._ibrion = 2

        self.set = """
LCHARG = .FALSE.
LWAVE  = .FALSE.
LELF   = .FALSE.
LVTOT  = .FALSE.
        """
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
        if self.in_type in ['unrelax', 'scf']:
            with open(fname, 'w') as fid:
                fid.write("""
PREC   = A
IBRION =-1
ISMEAR =-5
ISIF   = 2

ENCUT  = {:d}
EDIFF  = {:1.0e}
NPAR   = {:d}
{}
              """.format(self._encut,
                         self._accuracy,
                         self._npar,
                         self.set))
                fid.close()

        elif self.in_type in ['isif4']:
            with open(fname, 'w') as fid:
                fid.write("""
PREC   = A
NSW    = 100
IBRION = 2
ISIF   = 4

ISMEAR = 1
ENCUT  = 350.0
EDIFF  = 1e-05
NPAR   = 4
{}
                    """.format(self.set))

        elif self.in_type == 'dftrelax' or self.in_type == 'relax':
            with open(fname, 'w') as fid:
                fid.write("""
NSW    = 100
IBRION = 2
ISIF   = 2
ISMEAR = 1
NPAR   = 4
{}
""" .format(self.set))
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
                         """.format(self.nsw,
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
NELM   = 400
NBANDS = 118
{}
""".format(self.set))
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
