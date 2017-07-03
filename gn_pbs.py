#!/usr/bin/env python
# encoding: utf-8
#
###################################################################
#
# File Name : ./gn_pbs.py
#
###################################################################
#
# Purpose :     generate the  pbs file
#
# Creation Date :
# Last Modified : Sat Apr  1 23:15:50 2017
# Created By    : Chaoming Yang
#
###################################################################
import os


class gn_pbs(object):

    def __init__(self,
                 in_type="va",
                 job_title=None,
                 nnodes=2,
                 walltime=80,
                 exe_name=None):

        self.in_type = in_type
        self._nnodes = nnodes
        self._wall_time = walltime
        self._ppn = 12
        return

    def set_nnodes(self, nnodes=2):
        self._nnodes = nnodes
        return

    def get_nnodes(self):
        return self._nnodes

    def walltime(self):
        return self._wall_time

    def set_wall_time(self, wtime=12):
        self._wall_time = wtime
        return

    def get_wall_time(self):
        return self._wall_time

    def ppn(self):
        return self._ppn

    def set_ppn(self, ppn):
        self.ppn = ppn
        return

    def set_main_job(self, jobname):
        self.exe = jobname
        return

    def set_job_title(self, title):
        self.job_title = title
        return

    def set_pbs_type(self, in_type):
        self.in_type = in_type
        return

    def get_job_name(self):
        # try join
        list = ["A", "B", "C"]
        jobname = '-'.join(list)
        print jobname
        return

    def copy_inputs(self, dirname, *args):
        for filename in args:
            os.system("cp {} {}".format(filename, dirname))
        return

    def move_inputs(self, dirname, *args):
        for filename in args:
            os.system("mv {} {}".format(filename, dirname))
        return

    def write_pbs(self, od=None):
        if self.in_type == "va":
            print "assign_outfile"
            outfile = "va.pbs"
        if od is True:
            flux_type = "fluxod"
        else:
            flux_type = "flux"

        with open(outfile, 'w') as fid:
            fid.write("""#!/bin/sh
####  PBS preamble

#PBS -N %s
#PBS -M chaomy@umich.edu
#PBS -m e

#PBS -l nodes=%d:ppn=%d:ib,pmem=2gb,walltime=%d:00:00
#PBS -j oe
#PBS -V

#PBS -A qiliang_%s
#PBS -q %s
#PBS -l qos=flux

####  End PBS preamble
cd $PBS_O_WORKDIR
""" % (self.job_title, self._nnodes, self._ppn,
                self._wall_time, flux_type, flux_type))
            if type(self.exe) is list:
                fid.writelines(self.exe)
            else:
                fid.write("""{}""".format(self.exe))
            fid.close()
        return


if __name__ == '__main__':
    from optparse import OptionParser
    usage = "usage:%prog [options] arg1 [options] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option(
        "-t",
        "--type",
        action="store",
        type="string",
        dest="type",
        default="lattice")

    (options, args) = parser.parse_args()
