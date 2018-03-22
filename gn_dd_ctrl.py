#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-14 21:41:11

# Be #
from collections import OrderedDict
from .gn_dd_ctrl_hcp import *
from .gn_dd_ctrl_bcc import *


class gn_dd_ctrl(object):

    def __init__(self):
        self.ctrlfile = '{}.ctrl'.format(self.job)
        self.divisionline = '\n# ---------------------------------------------- #\n\n'
        self.fmat = "{:30} = {:18}\n"
        self.strmat = "{:30} = {:>18}\n"
        self.space = "\n"

    def set_fname(self, jobid):
        self.job = jobid
        filesetting['dirname'] = '{}_results'.format(self.job)
        filesetting['inclusionFile'] = '{}.dat'.format(self.job)
        self.ctrlfile = '{}.ctrl'.format(self.job)
        self.precfile = '{}.dat'.format(self.job)
        self.datafile = '{}.data'.format(self.job)
        print(self.ctrlfile)

    def write_mobility(self, fid):
        mobs = hcpmobs['Mg']
        types = ['a', 'cpa', 'c', 'o']
        # mobs = bccmobs['Fe']
        # types = ['g']
        for tkey in types:
            for key in list(mobs[tkey].keys()):
                fid.write(self.fmat.format(key, mobs[tkey][key]))
            fid.write(self.space)
        return fid

    def write_loading(self, fid, ltype='stress'):
        fid.write(self.divisionline)
        if ltype in ['stress']:
            key = 'appliedStress'
            stress = self.cal_stress(scale=-2.0e7)   # 10 to 20 Mpa
            fid.write('{:30} = '.format(key))
            fid.write('[{} {} {} {} {} {}]\n'.format(
                stress[0], stress[1], stress[2],
                stress[3], stress[4], stress[5]
            ))
        elif ltype in ['strain']:
            key = 'edotdir'
            fid.write('{:30} = '.format(key))
            fid.write(' [{:>5}{:>5}{:>5}]\n'.format(
                loading[key][0], loading[key][1], loading[key][2]
            ))
            key = 'eRate'
            fid.write(self.fmat.format(key, loading[key]))
        return fid

    def write_control(self, fid):
        settings = [domsetting, filesetting, matsetting['Mg'],
                    discretization, integration, nfmm, enforceglide]
        for ctrl in settings:
            fid.write(self.divisionline)
            for key in ctrl:
                if type(ctrl[key]) is str:
                    fid.write(self.strmat.format(key, ctrl[key]))
                else:
                    fid.write(self.fmat.format(key, ctrl[key]))
        return fid

    def write_outsetting(self, fid):
        ptype = ['flux', 'cn', 'prop', 'visit']
        fid.write(self.divisionline)
        for pkey in ptype:
            for key in outsetting[pkey]:
                fid.write(self.fmat.format(key, outsetting[pkey][key]))
            fid.write('\n')
        return fid

    def write_ctrl_file(self, ltype='stress'):
        fid = open(self.ctrlfile, 'w')
        fid = self.write_mobility(fid)
        fid = self.write_loading(fid, ltype=ltype)
        fid = self.write_control(fid)
        fid = self.write_outsetting(fid)
        fid.close()


if __name__ == '__main__':
    drv = gn_dd_ctrl()
    drv.cal_stress()
