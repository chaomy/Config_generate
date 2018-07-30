#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-29 16:50:26
# @Last Modified by:   chaomy
# @Last Modified time: 2018-07-01 02:19:09


class gn_md_input_thermo(object):

    def __init__(self, **kwargs):
        self.basics = kwargs

    def _write_equilibrium(self, **kwargs):
        basics = self.basics
        defaults = dict([('tstart', 0.0),
                         ('tend', 100),
                         ('tdamp', 0.1),
                         ('pstart', 0.0),
                         ('pend', 1.0),
                         ('pdamp', 1.0),
                         ('run1', 100000)])

        if kwargs is not None:
            for key, value in kwargs.items():
                if key in defaults:
                    defaults[key] = value

        with open(basics['innpt'], 'w') as fid:
            fid.write("""units          metal
atom_style     atomic

read_data      ../%s

pair_style     meam/spline    
pair_coeff     * *    ../meam.lib.best    Nb 
mass           *       92.906

neigh_modify   every   5      delay   0   check   yes

fix            100    all     box/relax   iso   0.0   vmax   0.1
minimize       1e-20  1e-20   1000000   1000000

write_restart  ../init_restart
                    """ % (basics['read_data']))
            fid.close()

    def _write_run_npt(self, **kwargs):
        basics = self.basics
        defaults = dict([('tstart', 0.1),
                         ('tend', 100),
                         ('tdamp', 0.1),
                         ('pstart', 0.0),
                         ('pend', 1.0),
                         ('pdamp', 1.0),
                         ('run1', 100000)])

        if kwargs is not None:
            for key, value in kwargs.items():
                if key in defaults:
                    defaults[key] = value

        with open(basics['innpt'], 'w') as fid:
            fid.write("""units          metal
read_data      lmp_init.txt

pair_style     meam/spline    
pair_coeff     * *    ../meam.lib.best    Nb 
mass           *      92.906

thermo         500
thermo_style   custom   step   etotal  pe  temp  press   pxx  pyy  pzz  lx  ly  lz

velocity       all      create  %g       4928459
fix            1        all     npt      temp     %g  %g  50   iso  %g  %g  100 drag 1.0
#dump          1        all     custom   10000    dump.*  id  type  x y z
restart        500000   rst.* 
run            5000000
unfix          1 
                    """ % (2.0 * defaults['tend'], defaults['tend'], defaults['tend'],
                           defaults['pend'], defaults['pend']))
            fid.close()

        with open('in.rst', 'w') as fid:
            fid.write("""units          metal
read_restart   rst.500000 

pair_style     meam/spline    
pair_coeff     * *    ../meam.lib.best    Nb 
mass           *      92.906

thermo         500
thermo_style   custom   step   etotal  pe  temp  press   pxx  pyy  pzz  lx  ly  lz

fix            1        all     npt      temp     %g  %g  50   iso  %g  %g  100 drag 5.0
restart        500000   rst.* 
run            2500000
unfix          1 
                    """ % (defaults['tend'], defaults['tend'], defaults['pend'], defaults['pend']))
            fid.close()
