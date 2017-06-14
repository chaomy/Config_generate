#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : gn_md_input_perf.py
#
###################################################################
#
# Purpose :
#
# Creation Date :
# Last Modified : Tue Apr 18 12:45:39 2017
# Created By    : Chaoming Yang
#
###################################################################


class gn_md_input_cu3au(object):
    def __init__(self):

        return

    def _write_Cu3Au_restart(self, fname):
        with open(fname, 'w') as fid:
            fid.write("""# ---------- Initialize Simulation ---------------------
clear
units       metal
atom_style  atomic
boundary    p  p  p
atom_modify sort  0 0.0

# ---------- Create Atomistic Structure ------------
read_data   ./lmp_init.txt

# ---------- Define Interatomic Potential -------------------
pair_style   eam/alloy
pair_coeff   *  *  ../CuAg.eam.alloy.txt  Cu  Ag
neighbor     1.0    bin
neigh_modify every  10  delay 0 check yes

#delete_atoms overlap   0.3   all   all
# ---------- Run Minimization ---------------------
thermo       1000
thermo_style  custom  step  etotal  pe  temp  lx  ly  lz
dump       mcust   all  custom   300   custom/md.*.dump  id  type  mass x y z
dump_modify  mcust   pad    5

min_style cg
minimize  1e-12   1e-10   100000   100000
write_restart     restart/init_restart

#restart     5000    restart/restart.*
                    """)
            fid.close()
        return

    def _write_Cu3Au_perf(self, fname, lattice, size):
        with open(fname, 'w') as fid:
            fid.write("""# ---------- Initialize Simulation ---------------------
clear
units        metal
atom_style   atomic
boundary     p   p   p
atom_modify  sort  0 0.0

# ---------- Create Atomistic Structure ------------
variable  alat   equal    %f 
variable  xdim   equal    %d*${alat}*sqrt(2)/2.
variable  yshif  equal    10*${alat}*sqrt(2)/2.
variable  ydim   equal    (${yshif}+%d*${alat}*sqrt(2)/2.)
variable  nydim  equal    (${yshif}-3*${alat}*sqrt(2)/2.)
variable  uydim  equal    (${yshif}+60*${alat}*sqrt(2)/2.)
variable  lydim  equal    (${yshif}-${yshif})
variable  uinter equal    (0.05+${yshif})
variable  linter equal    (-0.05+${yshif})
variable  fixup  equal    ${linter}+20
variable  zdim   equal    %d*${alat}

#----------- Params Fcc ----------------------------
region       whole block    0.000   ${xdim}  ${lydim}  ${uydim}   0.000   ${zdim}
create_box   4  whole

#----------- Params L12 ----------------------------
region   rCu3Au   block    INF   INF   ${nydim}  ${uinter}  INF   INF   units  box
lattice  custom  ${alat}  &
a1   1.0  0.0  0.0  &
a2   0.0  1.0  0.0  &
a3   0.0  0.0  1.0  &
basis  0.0   0.0   0.0   &
basis  0.0   0.5   0.5   &
basis  0.5   0.5   0.0   &
basis  0.5   0.0   0.5   &
orient  x   -1   1   0  &
orient  y    1   1   0  &
orient  z    0   0  -1

create_atoms  1  region  rCu3Au &
basis  1    1    &
basis  2    1    &
basis  3    1    &
basis  4    2

region   rPureCu block   INF   INF    ${linter}   ${ydim}   INF  INF units  box 
lattice   custom  ${alat} &
a1  1.0  0.0  0.0  &
a2  0.0  1.0  0.0  &
a3  0.0  0.0  1.0  &
basis  0.0   0.0   0.0   &
basis  0.0   0.5   0.5   &
basis  0.5   0.5   0.0   &
basis  0.5   0.0   0.5   &
orient  x   -1   1    0  &
orient  y    1   1    0  &
orient  z    0   0   -1

create_atoms  3  region  rPureCu &
basis  1    3    &
basis  2    3    &
basis  3    3    &
basis  4    3

region  fix    block  INF   INF   ${fixup}   ${ydim}   INF   INF    units   box
group   fix    region  fix
set     group  fix    type  4

# ---------- Define Interatomic Potential -------------------
pair_style  eam/alloy
pair_coeff  *  *  ../CuAg.eam.alloy.txt  Cu  Ag   Cu   Cu
neighbor   1.0  bin
neigh_modify  every  10  delay 0 check yes

delete_atoms   overlap   0.3   all   all

dump     mcust   all  custom  500   perf.dump  id  type  mass x y z
thermo   1000
thermo_style  custom  step  pe  temp  lx  ly  lz

run   1
                    """ % (lattice, size.xdim, size.ydim, size.zdim))
            fid.close()
        return
