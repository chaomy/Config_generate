#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name :
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified :
# Created By    : Chaoming Yang
#
###################################################################


class gn_md_input_thermo(object):
    def __init__(self, **kwargs):
        self.basics = kwargs
        return

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
            for key, value in kwargs.iteritems():
                if key in defaults:
                    defaults[key] = value

        with open(basics['innpt'], 'w') as fid:
            fid.write("""clear
units          metal
dimension      3
boundary       p p p
atom_style     atomic

read_data      ../%s

pair_style     %s
pair_coeff     * *    ../%s   %s
neighbor       2.0     bin
neigh_modify   every   5    delay   0   check   yes

fix            100    all     box/relax   iso   0.0   vmax   0.1
minimize       1e-20  1e-20   1000000   1000000

write_restart  ../init_restart
                    """ % (basics['read_data'],
                           basics['pair_style'],
                           basics['pair_coeff'],
                           basics['element']))
            fid.close()
        return

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
            for key, value in kwargs.iteritems():
                if key in defaults:
                    defaults[key] = value

        with open(basics['innpt'], 'w') as fid:
            fid.write("""clear
units          metal
dimension      3
boundary       p p p
atom_style     atomic

read_restart   ../%s
reset_timestep  0

pair_style     %s
pair_coeff     * *     ../%s   %s
neighbor       2.0     bin
neigh_modify   every   5    delay   0   check   yes

timestep       0.001

thermo         500
thermo_style   custom   step   etotal   pe   temp   press   lx  ly  lz
#dump          mcust   all  custom   10000    out/md.*.dump  id  type  mass x y z

fix            1       all     npt     temp     %g    %g      0.1    iso   %g    %g    1.0

run            %d
unfix          1

fix            2       all     npt     temp     %g     %g     0.1    iso   %g    %g    1.0
run            300000
unfix          2
write_restart  ../init_restart
                    """ % (basics['read_restart'],
                           basics['pair_style'],
                           basics['pair_coeff'],
                           basics['element'],
                           defaults['tstart'],
                           defaults['tend'],
                           defaults['pstart'],
                           defaults['pend'],
                           defaults['run1'],
                           defaults['tend'],
                           defaults['tend'],
                           defaults['pend'],
                           defaults['pend']))
            fid.close()
        return
