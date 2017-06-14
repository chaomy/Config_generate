#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : ../Generate_config/gn_md_input_fcc_lattice.py
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


class gn_md_input_fcc_lattice(object):
    def __init__(self):

        return

    def _write_fcc_lattice(self,
                           fname,
                           lattice,
                           element,
                           pottype,
                           potname):
        with open(fname, 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
variable     alat   equal   %f
region       whole block    0.000   ${alat}  0  ${alat}   0.0   ${alat}
create_box   1   whole

lattice        fcc   ${alat}
create_atoms   1  region  whole
 
# --------------------- FORCE FIELDS ---------------------
pair_style       %s
pair_coeff       *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    5    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute       peratom     all     pe/atom
compute       eatoms      all     reduce   sum     c_peratom
thermo            1000

thermo_style      custom  step  atoms etotal  pe   c_eatoms  lx   ly   lz

variable     Ef     equal   "c_eatoms"

#--------------------- MINIMIZE -------------------------------
fix          1    all    box/relax  x  0.0  y  0.0  z  0.0  couple none  vmax 0.001
min_style    cg
minimize     1e-20      1e-20     100000     100000
run          1
unfix        1

variable lengthx  equal "lx"
variable lengthy  equal "ly"
variable lengthz  equal "lz"
variable EperAtom equal "etotal/atoms"

print "Lattice_constantx = ${lengthx}"
print "Lattice_constanty = ${lengthy}"
print "Lattice_constantz = ${lengthz}"
print "energy per atom = ${EperAtom} eV/atom"
                    """ % (lattice, pottype, potname, element))
        return
