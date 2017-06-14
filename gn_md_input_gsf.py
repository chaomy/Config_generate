#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : gn_md_input_gsf.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified : Mon Apr 24 12:56:42 2017
# Created By    : Chaoming Yang
#
###################################################################


class gn_md_input_gsf(object):
    def __init__(self):
        return

    def _write_gsf_minimize(self,
                            config_file,
                            potential_type,
                            potential_file,
                            element_name,
                            tag):
        if tag == 'relaxed':
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
                        """ % (config_file,
                               potential_type,
                               potential_file,
                               element_name))
        elif tag == 'unrelaxed':
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
                        """ % (config_file,
                               potential_type,
                               potential_file,
                               element_name))
        return
