#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2018-05-28 13:15:01

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
        self._mem = 2 

    def set_nnodes(self, nnodes=2):
        self._nnodes = nnodes

    def get_nnodes(self):
        return self._nnodes

    def walltime(self):
        return self._wall_time

    def set_wall_time(self, wtime=12):
        self._wall_time = wtime

    def get_wall_time(self):
        return self._wall_time

    def ppn(self):
        return self._ppn

    def set_ppn(self, ppn):
        self._ppn = ppn
        
    def set_mem(self, mem):
        self._mem = mem

    def set_main_job(self, jobname):
        self.exe = jobname

    def set_job_title(self, title):
        self.job_title = title

    def set_pbs_type(self, in_type):
        self.in_type = in_type

    def get_job_name(self):
        list = ["A", "B", "C"]
        jobname = '-'.join(list)
        print(jobname)

    def copy_inputs(self, dirname, *args):
        for filename in args:
            os.system("cp {} {}".format(filename, dirname))

    def move_inputs(self, dirname, *args):
        for filename in args:
            os.system("mv {} {}".format(filename, dirname))

    def write_pbs(self, od=False):
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

#PBS -l nodes=%d:ppn=%d:ib,pmem=%dgb,walltime=%d:00:00
#PBS -j oe
#PBS -V

#PBS -A qiliang_%s
#PBS -q %s
#PBS -l qos=flux

####  End PBS preamble
cd $PBS_O_WORKDIR
""" % (self.job_title, self._nnodes, self._ppn, self._mem,
                self._wall_time, flux_type, flux_type))
            if type(self.exe) is list:
                fid.writelines(self.exe)
            else:
                fid.write("""{}""".format(self.exe))
            fid.close()


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
