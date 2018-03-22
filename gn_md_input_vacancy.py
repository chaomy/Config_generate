#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2018-03-20 08:49:15
# @Last Modified by:   chaomy
# @Last Modified time: 2018-03-20 08:49:28


class gn_md_input_vacancy(object):

    def _write_input_vacancy(self,
                             config_file,
                             potential_type,
                             potential_file,
                             element_name):
        with open("in.minimize", 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data     %s
# --------------------- FORCE FIELDS ---------------------
pair_style      %s
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    5    delay   0
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom
thermo            1000
thermo_style      custom   step   pe   c_eatoms   etotal

variable     Ef     equal   "c_eatoms"

dump     1    all    custom    10000   bcc.init.*  id type mass  x y z
#--------------------- MINIMIZE -------------------------------
min_style    cg
minimize     1e-16      1e-16     100000     100000
                    """ % (config_file,
                           potential_type,
                           potential_file,
                           element_name))

    def _write_neb_init(self,
                        config_file,
                        potential_type,
                        potential_file,
                        element_name):
        with open("in.init", 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
atom_modify	  map   array   sort   0  0.0
# --------------------- FORCE FIELDS ---------------------
read_data     %s
# --------------------- FORCE FIELDS ---------------------
pair_style      %s
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    5    delay   0
# --------------------- SETTINGS ---------------------
thermo            1000
thermo_style      custom   step   pe   etotal

dump     1    all    custom    10000   bcc.init.*  id type mass  x y z
dump_modify   1       pad   5

#--------------------- MINIMIZE -------------------------------
min_style       cg
minimize        1e-15      1e-15     100000     100000
write_restart   init_restart
                    """ % (config_file,
                           potential_type,
                           potential_file,
                           element_name))

    def _write_neb_final(self,
                         config_file,
                         potential_type,
                         potential_file,
                         element_name):
        with open("in.final", 'w') as fid:
            fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units		metal
dimension	3
boundary	p   p   p
atom_style	atomic
atom_modify	  map   array   sort   0  0.0

# --------------------- ATOM DEFINITION ---------------------#
read_data  %s

pair_style  %s 
pair_coeff  *  *  %s  %s 
neighbor    2.0     bin
neigh_modify    every   5   delay   0   check yes

#--------------------set neb group--------------
thermo    100

dump      1   all   custom   10000   bcc.final.*  id type mass  x y z
dump_modify   1       pad   5

min_style        cg
minimize         1e-15  1e-15  100000  100000
                    """ % (config_file,
                           potential_type,
                           potential_file,
                           element_name))

    def _write_neb_run(self,
                       config_file,
                       potential_type,
                       potential_file,
                       element_name):
        with open("in.neb", 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units		metal
dimension	3
boundary	p   p   p
atom_style	atomic
atom_modify	  map   array   sort   0  0.0
variable      Ulabel uloop   24

# --------------------- ATOM DEFINITION ---------------------#
read_restart  %s

pair_style  %s
pair_coeff  *  *  %s  %s 
neighbor    2.0     bin
neigh_modify    every   5   delay   0   check yes

#--------------------set neb group--------------
thermo    100

fix     1     all    neb        3.0

dump    1     all    custom   10000   dump.all.${Ulabel}  id type mass  x y z
dump_modify  1  pad  5

min_style     fire
neb     0.0   0.02    100000    100000    200     final   final.coord
""" % (config_file, potential_type, potential_file, element_name))
