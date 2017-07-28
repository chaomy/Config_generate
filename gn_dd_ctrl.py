#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-07-27 20:31:25


hcp_a_mobility = {
    "HCP_A_Basal_EdgeDrag": 1.18434387898773e-04,
    "HCP_A_Basal_ScrewDrag": 1.65421659553579e-04,
    "HCP_A_Prismatic_EdgeDrag": 5.69513923115224e-05,
    "HCP_A_Prismatic_ScrewDrag": 7.15096201802962e-05,
    "HCP_A_1stPyramidal_EdgeDrag": 7.28268834977848e-05,
    "HCP_A_1stPyramidal_ScrewDrag": 1.05154866364405e-04,
    "HCP_A_2ndPyramidal_EdgeDrag": 7.28268834977848e-05,
    "HCP_A_2ndPyramidal_ScrewDrag": 1.05154866364405e-04,
}
hcp_apc_mobility = {
    # Assumed 20 times basal edge
    "HCP_CpA_Prismatic_EdgeDrag": 2.36868775797546e-03,
    "HCP_CpA_Prismatic_ScrewDrag": 2.36868775797546e-03,
    "HCP_CpA_1stPyramidal_EdgeDrag": 2.36868775797546e-03,
    "HCP_CpA_1stPyramidal_ScrewDrag": 2.36868775797546e-03,
    "HCP_CpA_2ndPyramidal_EdgeDrag": 2.36868775797546e-03,
    "HCP_CpA_2ndPyramidal_ScrewDrag": 2.36868775797546e-03,
}
hcp_c_mobility = {
    "HCP_C_Prismatic_EdgeDrag": 2.36868775797546e-03,
    "HCP_C_Prismatic_ScrewDrag": 2.36868775797546e-03,
    "HCP_Sessile_ScrewDrag": 23.6868,  # 10000 * pyramidal
    "HCP_Sessile_EdgeDrag": 23.6868,  # 10000 * pyramidal
    "HCP_LineDrag": 5.e-1
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

domsetting = {
    'numXdoms': 1,
    'numYdoms': 1,
    'numZdoms': 1,
    'numXcells': 4,
    'numYcells': 4,
    'numZcells': 4
}

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

outsetting = {
    'fluxfile': 1,
    'fluxfreq': 200,
    'fluxdt': -1.0,
    'fluxtime': 0.0,
    'fluxcounter': 0,
    'savecn': 1,
    'savecnfreq': 200,
    'savecndt': -1.0,
    'savecntime': 0.0,
    'savecncounter': 0,
    'saveprop': 1,
    'savepropfreq': 200,
    'savepropdt': -1.0,
    'saveproptime': 0.0,
    'writeVisit': 1,
    'writeVisitFreq': 200,
    'writeVisitCounter': 1,
    'writeVisitSegments': 1,
    'writeVisitNodes': 1,
    'writeVisitSegmentsAsText': 1,
    'writeVisitNodesAsText': 1
}


class gn_dd_ctrl(object):

    def __init__(self):
        self.ctrlfile = '{}.ctrl'.format(self.job)
        self.divisionline = '\n# ---------------------------------------------- #\n\n'
        self.fmat = "{:30} = {:18}\n"
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
        for key in hcp_a_mobility:
            fid.write(self.fmat.format(key, hcp_a_mobility[key]))
        fid.write(self.space)
        for key in hcp_apc_mobility:
            fid.write(self.fmat.format(key, hcp_apc_mobility[key]))
        fid.write(self.space)
        for key in hcp_c_mobility:
            fid.write(self.fmat.format(key, hcp_c_mobility[key]))
        fid.write(self.space)
        return fid

    def write_loading(self, fid, opt='stress'):
        fid.write(self.divisionline)
        if opt in ['stress']:
            key = 'appliedStress'
            stress = self.cal_stress(scale=1e8)
            fid.write('{:30} = '.format(key))
            fid.write('[{} {} {} {} {} {}]\n'.format(
                stress[0], stress[1], stress[2],
                stress[3], stress[4], stress[5]
            ))
        elif opt in ['strain']:
            key = 'edotdir'
            fid.write('{} = '.format(key))
            fid.write('[{} {} {}] \n'.format(
                loading[key][0], loading[key][1], loading[key][2]
            ))
            key = 'eRate'
            fid.write('{} = {}'.format(key, loading[key]))
        return fid

    def write_control(self, fid):
        settings = [domsetting, filesetting, matsetting,
                    discretization, integration, fmm, outsetting]
        for ctrl in settings:
            fid.write(self.divisionline)
            for key in ctrl:
                fid.write(self.fmat.format(key, ctrl[key]))
        return fid

    def write_ctrl_file(self):
        fid = open(self.ctrlfile, 'w')
        fid = self.write_hcp_mobility(fid)
        fid = self.write_loading(fid)
        fid = self.write_control(fid)
        fid.close()
        return


if __name__ == '__main__':
    drv = gn_dd_ctrl()
    drv.cal_stress()
