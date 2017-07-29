#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-28 16:41:20


# Be #
from collections import OrderedDict

hcpmobs = {
    'Be':
    {
        'a': OrderedDict(
            [("HCP_A_Basal_EdgeDrag", 1.18434387898773e-04),
             ("HCP_A_Basal_ScrewDrag", 1.65421659553579e-04),
             ("HCP_A_Prismatic_EdgeDrag", 5.69513923115224e-05),
             ("HCP_A_Prismatic_ScrewDrag", 7.15096201802962e-05),
             ("HCP_A_1stPyramidal_EdgeDrag", 7.28268834977848e-05),
             ("HCP_A_1stPyramidal_ScrewDrag", 1.05154866364405e-04),
             ("HCP_A_2ndPyramidal_EdgeDrag", 7.28268834977848e-05),
             ("HCP_A_2ndPyramidal_ScrewDrag", 1.05154866364405e-04)
             ]),  # Assumed 20 times basal edge
        'cpa': OrderedDict([
            ("HCP_CpA_Prismatic_EdgeDrag", 2.36868775797546e-03),
            ("HCP_CpA_Prismatic_ScrewDrag", 2.36868775797546e-03),

            ("HCP_CpA_1stPyramidal_EdgeDrag", 2.36868775797546e-03),
            ("HCP_CpA_1stPyramidal_ScrewDrag", 2.36868775797546e-03),

            ("HCP_CpA_2ndPyramidal_EdgeDrag", 2.36868775797546e-03),
            ("HCP_CpA_2ndPyramidal_ScrewDrag", 2.36868775797546e-03)
        ]),  # Assumed 20 times basal edge
        'c': OrderedDict([
            ("HCP_C_Prismatic_EdgeDrag", 2.36868775797546e-03),
            ("HCP_C_Prismatic_ScrewDrag", 2.36868775797546e-03)
        ]),
        'o': OrderedDict([
            ("HCP_Sessile_ScrewDrag", 23.6868),  # 10000 * pyramidal
            ("HCP_Sessile_EdgeDrag", 23.6868),  # 10000 * pyramidal
            ("HCP_LineDrag", 5.e-6)
        ])},
    'Mg':
    {
        'a': OrderedDict([
            ("HCP_A_Basal_EdgeDrag", 4.7e-6),
            ("HCP_A_Basal_ScrewDrag", 1.3e-5),

            ("HCP_A_Prismatic_EdgeDrag", 7.7e-6),
            ("HCP_A_Prismatic_ScrewDrag", 3.7e-5),

            ("HCP_A_1stPyramidal_EdgeDrag", 8.0e-5),
            ("HCP_A_1stPyramidal_ScrewDrag", 8.0e-5),

            ("HCP_A_2ndPyramidal_EdgeDrag", 8.0e-5),
            ("HCP_A_2ndPyramidal_ScrewDrag", 8.0e-5)]),
        # Assumed 20 times basal edge
        'cpa': OrderedDict([
            ("HCP_CpA_Prismatic_EdgeDrag", 9.4e-05),
            ("HCP_CpA_Prismatic_ScrewDrag", 9.4e-05),

            ("HCP_CpA_1stPyramidal_EdgeDrag", 9.4e-05),
            ("HCP_CpA_1stPyramidal_ScrewDrag", 9.4e-05),

            ("HCP_CpA_2ndPyramidal_EdgeDrag", 9.4e-05),
            ("HCP_CpA_2ndPyramidal_ScrewDrag", 9.4e-05)]),
        # Assumed 20 times basal edge
        'c': OrderedDict([
            ("HCP_C_Prismatic_EdgeDrag", 9.4e-05),
            ("HCP_C_Prismatic_ScrewDrag", 9.4e-05)]),
        'o': OrderedDict([
            ("HCP_Sessile_ScrewDrag", 0.94),  # 10000 * pyramidal
            ("HCP_Sessile_EdgeDrag", 0.94),  # 10000 * pyramidal
            ("HCP_LineDrag", 5.e-6)
        ])
    }
}

loading = {
    'loadType': 0,
    'edotdir': [1.0, 0.0, 0.0],
    'eRate': 1e5,
    'appliedStress': None
}

matsetting = {
    'rc': 3.0,
    'pois': 0.29,
    'shearModulus': 45e9,
    'burgMag': 3.2094e-10,
    'cOVERa': 1.6236,
    'mobilityLaw': 'HCP_linear_Eshelby'
}

filesetting = {
    'dirname': 'HCP_results',
    'inclusionFile': 'prect.dat'
}

domsetting = OrderedDict([
    ('numXdoms', 1),
    ('numYdoms', 1),
    ('numZdoms', 1),
    ('numXcells', 4),
    ('numYcells', 4),
    ('numZcells', 4)])


discretization = {
    'remeshRule': 2,
    'maxSeg': 200.,
    'minSeg': 50.
}

integration = {
    'maxstep': 1000000000,
    'timestepIntegrator': "trapezoid",
    'nextDT': 1.0e-12,
    'rTol': 1.5,
    'rann': 1.0,
    'collisionMethod': 1
}

fmm = {
    'fmEnabled': 0,
    'fmMPOrder': 2,
    'fmTaylorOrder': 5,
    'fmCorrectionTbl': "fm-ctab.hcp.m2.t5.dat"
}

freq = 300
outsetting = {
    'flux': {
        'fluxfile': 1,
        'fluxfreq': freq,
        'fluxdt': -1.0,
        'fluxtime': 0.0,
        'fluxcounter': 0,
    },
    'cn': {
        'savecn': 1,
        'savecnfreq': freq,
        'savecndt': -1.0,
        'savecntime': 0.0,
        'savecncounter': 0,
    },
    'prop': {
        'saveprop': 1,
        'savepropfreq': freq,
        'savepropdt': -1.0,
        'saveproptime': 0.0,
    },
    'visit': {
        'writeVisit': 1,
        'writeVisitFreq': freq,
        'writeVisitCounter': 1,
        'writeVisitSegments': 1,
        'writeVisitNodes': 1,
        'writeVisitSegmentsAsText': 1,
        'writeVisitNodesAsText': 1
    }
}


class gn_dd_ctrl(object):

    def __init__(self):
        self.ctrlfile = '{}.ctrl'.format(self.job)
        self.divisionline = '\n# ---------------------------------------------- #\n\n'
        self.fmat = "{:30} = {:18}\n"
        self.strmat = "{:30} = {:>18}\n"
        self.space = "\n"
        return

    def set_fname(self, jobid):
        self.job = jobid
        filesetting['dirname'] = '{}_results'.format(self.job)
        filesetting['inclusionFile'] = '{}.dat'.format(self.job)
        self.ctrlfile = '{}.ctrl'.format(self.job)
        self.precfile = '{}.dat'.format(self.job)
        self.datafile = '{}.data'.format(self.job)
        print self.ctrlfile
        return

    def write_hcp_mobility(self, fid):
        mobs = hcpmobs['Mg']
        types = ['a', 'cpa', 'c', 'o']
        for tkey in types:
            for key in mobs[tkey].keys():
                fid.write(self.fmat.format(key, mobs[tkey][key]))
            fid.write(self.space)
        return fid

    def write_loading(self, fid, ltype='stress'):
        fid.write(self.divisionline)
        if ltype in ['stress']:
            key = 'appliedStress'
            stress = self.cal_stress(scale=1e8)
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
        settings = [domsetting, filesetting, matsetting,
                    discretization, integration, fmm]
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
        fid = self.write_hcp_mobility(fid)
        fid = self.write_loading(fid, ltype=ltype)
        fid = self.write_control(fid)
        fid = self.write_outsetting(fid)
        fid.close()
        return


if __name__ == '__main__':
    drv = gn_dd_ctrl()
    drv.cal_stress()
