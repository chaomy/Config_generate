# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-12-09 22:26:05
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-14 22:15:50

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
            ("HCP_LineDrag", 1e-5)
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
            ("HCP_LineDrag", 1e-5),
            ("MobEshelbyResist", 1e5)
        ])
    }
}

loading = {
    'loadType': 0,
    'edotdir': [1.0, 0.0, 0.0],
    'eRate': 1e5,
    'appliedStress': None
}

matsetting = {'Mg': {'rc': 3.0,
                     'pois': 0.29,
                     'shearModulus': 17e9,
                     'burgMag': 3.2094e-10,
                     'cOVERa': 1.6236,
                     'mobilityLaw': 'HCP_linear_Eshelby'},
              'Fe': {'rc': 3.0,
                     'pois': 0.291,
                     'shearModulus': 94e9,
                     'burgMag': 2.4727727e-10,
                     'mobilityLaw': 'BCC_Linear'}
              }

filesetting = {
    'dirname': 'results',
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

enforceglide = {
    'enforceGlidePlanes': 1,
    'enableCrossSlip': -1,
    'elasticinteraction': 1
}

nfmm = {
    'fmEnabled': 0,
    'Rijmfile': "Rijm.cube.out",
    'RijmPBCfile': "RijmPBC.cube.out"
}


fmm = {
    'fmEnabled': 1,
    'fmMPOrder': 2,
    'fmTaylorOrder': 5,
    'fmCorrectionTbl': "fm-ctab.hcp.m2.t5.dat"
}

freq = 5000
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
