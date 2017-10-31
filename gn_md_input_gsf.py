#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   chaomy
# @Last Modified time: 2017-09-08 13:14:25


class gn_md_input_gsf(object):
    def _write_gsf_minimize(self, config_file,
                            potential_type,
                            potential_file,
                            element_name, tag):
        if tag in ['relaxed']:
            with open("in.md_gsf", 'w') as fid:
                fid.write("""
# --------------------- INITIALIZAITION ---------------------
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
    neigh_modify    every    1    delay   0
# --------------------- SETTINGS ---------------------
    thermo          1000
    thermo_style    custom   step    pe    etotal
    dump            1        all    custom   10000   bcc.init.*  id  type  mass  x y z
#--------------------- MINIMIZE -------------------------------
    fix          1     all  setforce    0    0    NULL
    min_style    cg
    minimize     1e-15     1e-15     100000     100000
    unfix        1
    variable  Etol  equal "etotal"
    print  "${Etol}"  file  out.txt
                        """ % (config_file,
                               potential_type,
                               potential_file,
                               element_name))
        elif tag in ['unrelaxed']:
            with open("in.md_gsf", 'w') as fid:
                fid.write("""
# --------------------- INITIALIZAITION ---------------------
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
    neigh_modify    every    1    delay   0
# --------------------- SETTINGS ---------------------
    thermo          1000
    thermo_style    custom   step    pe    etotal

#--------------------- MINIMIZE -------------------------------
    fix          1     all  setforce    0    0    0
    min_style    cg
    minimize     1e-18      1e-18     100000     100000
    unfix        1
    variable  Etol  equal "etotal"
    print  "${Etol}"  file  out.txt
                        """ % (config_file,
                               potential_type,
                               potential_file,
                               element_name))
        return
