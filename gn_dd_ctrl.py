#!/usr/local/bin/python
# encoding: utf-8


import os


class gn_dd_ctrl(object):

    def __init__(self):

        return

    def write_ctrl(self, screw_mob, edge_mob):
        with open("glide.ctrl", 'w') as fid:
            fid.write("""dirname = "."
#
numXdoms = 2
numYdoms = 2
numZdoms = 2
#
numXcells = 4
numYcells = 4
numZcells = 4
#
#  Load balancing setup
#
DLBfreq = 1
decompType = 2
#
#  Discretization controls #
remeshRule = 2
maxSeg = 5.000000e+02
enforceGlidePlanes = 1
enableCrossSlip = 1
#
mobilityLaw = "BCC_0"
#
maxstep = 10000000
timestepIntegrator = "trapezoid"
nextDT = 1.0e-08
#
#  FMM setup
#
fmEnabled       = 1
fmMPOrder       = 2
fmTaylorOrder   = 5
fmCorrectionTbl = "W_miu160_niu2p79"
#
#  Loading conditions
#
loadType = 1
eRate = 1.0e4  # unit 1/s
rc = 6

#  mobility parameters
pois = 2.796143e-01
shearModulus = 1.6e11
YoungModulus = 5.2e11
MobScrew = %e    # default 10
MobEdge =  %e    # default 10
MobClimb = 3.0e-02 # default 1.0e-02

#
savetimers = 0
savecn = 1
savecnfreq = 20000
savecncounter = 0

saveprop = 1
savepropfreq = 5000

writeVisit =       1
writeVisitFreq =   20000
writeVisitCounter =    0
writeVisitSegments =   1
writeVisitNodes =      1
                    """ % (screw_mob, edge_mob))
        return
