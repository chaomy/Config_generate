#!/usr/bin/env python
# encoding: utf-8
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-11-01 01:50:49
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-01 01:51:37


class gn_md_input_hcp_lattice(object):
    def _write_hcp_lattice(self,
                           fname,
                           lattice_a,
                           element,
                           pottype, potname):
        with open(fname, 'w') as fid:
            fid.write("""
# ---------- Initialize Simulation ---------------------
clear
units       metal
atom_style  atomic
boundary    p p p

# ---------- Create Atomistic Structure ---------------------
variable  alat  equal    %f
variable  clat  equal    sqrt(8./3.)*${alat} 
variable  xsize equal    1
variable  ysize equal    1
variable  zsize equal    1
variable  xdim  equal    ${alat}*${xsize}
variable  ydim  equal    sqrt(3)*${alat}*${ysize}
variable  zdim  equal    ${clat}*${zsize}
variable  sqa   equal    sqrt(3)
variable  b1    equal    1./3.
variable  b2    equal    5./6.
variable  c     equal    ${clat}/${alat}

region      whole  block   0   ${xdim}   0   ${ydim}    0    ${zdim}
create_box     1   whole

lattice   custom    ${alat} &
a1  1.0     0.0     0.0     &
a2  0.0     ${sqa}  0.0     &
a3  0.0     0.0     ${c}    &
basis   0.0   0.0   0.0   &
basis   0.5   0.5   0.0   &
basis   0.0   ${b1} 0.5   &
basis   0.5   ${b2} 0.5

create_atoms   1  box

# ---------- Define Interatomic Potential -------------------
pair_style     %s
pair_coeff     *  *   %s  %s

neighbor      1.0  bin
neigh_modify  every  10  delay 0 check yes

thermo   10
thermo_style  custom  step  atoms  etotal  pe  lx  ly  lz

# ---------- Run Minimization ---------------------
dump    mcfg  all  cfg   500  hcp.*.cfg  mass  type xs ys zs

variable fnlx  equal "lx"
variable fnly  equal "ly"
variable fnlz  equal "lz"
variable EperAtom equal "etotal/atoms"

min_style cg
fix      1   all box/relax  x  0.0  y  0.0  z  0.0  couple none
minimize 1e-18 1e-18 100000 100000

print "Lattice_constantx = ${fnlx}"
print "Lattice_constanty = ${fnly}"
print "Lattice_constantz = ${fnlz}"
print "energy per atom = ${EperAtom} eV/atom"
                    """ % (lattice_a,
                           pottype,
                           potname,
                           element))
            fid.close()
        return
