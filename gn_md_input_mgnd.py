#!/usr/bin/env python
# encoding: utf-8

###################################################################
#
# File Name : gn_md_input_mgnd.py
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


class gn_md_input_mgnd(object):
    def __init__(self):
        # 4 atom mg
        self.Eatom_mg_hcp = -6.1973714007776

        # a unit cell  16 (atoms)
        self.Ecoh_mg3_nd = -33.9878678899443
        return

    def _write_mg_nd_inter_phase(self, fname, param):
        with open(fname, 'w') as fid:
            fid.write("""# ---------- Initialize Simulation ---------------------
clear
units       metal
atom_style  atomic
boundary    p  s  p
atom_modify sort  0 0.0

# ---------- Params D03 ---------------------
variable  altD  equal    7.4662780330786
variable  xsize equal    5
variable  ysize equal    20
variable  zsize equal    4
variable  fixly equal    4

# ---------- Params Hcp ---------------------
variable  alat  equal    3.2019267694893
variable  clat  equal    5.1969105399
variable  xdim  equal    ${alat}*${xsize}
variable  ydim  equal    sqrt(3)*${alat}*${ysize}
variable  fydim equal    sqrt(3)*${alat}*(${ysize}-${fixly})
variable  nydim equal    -${ydim}
variable  nfydim equal   -${fydim}
variable  zdim  equal    ${clat}*${zsize}
variable  sqa   equal    sqrt(3)
variable  b1    equal    1./3.
variable  b2    equal    5./6.
variable  c     equal    ${clat}/${alat}

region    whole block    0.000   ${xdim}  ${nydim}  ${ydim}   0.000  ${zdim}
create_box   6   whole

region    rD03  block    INF  INF  ${nfydim}  0.05  INF  INF  units  box
lattice  custom  ${altD}  &
a1  1.0  0.0  0.0  &
a2  0.0  1.0  0.0  &
a3  0.0  0.0  1.0  &
basis  0.0   0.0   0.0   &
basis  0.0   0.5   0.5   &
basis  0.5   0.5   0.0   &
basis  0.5   0.0   0.5   &
basis  0.5   0.0   0.0   &
basis  0.0   0.5   0.0   &
basis  0.0   0.0   0.5   &
basis  0.5   0.5   0.5   &
basis  0.25  0.25  0.25  &
basis  0.75  0.25  0.25  &
basis  0.25  0.75  0.25  &
basis  0.25  0.25  0.75  &
basis  0.25  0.75  0.75  &
basis  0.75  0.25  0.75  &
basis  0.75  0.75  0.25  &
basis  0.75  0.75  0.75  &
orient  x   -1   1   -1  &
orient  y   -1   1    2  &
orient  z    1   1    0
create_atoms  2  region  rD03  &
basis    1    2   &
basis    2    2   &
basis    3    2   &
basis    4    2   &
basis    5    1   &
basis    6    1   &
basis    7    1   &
basis    8    1   &
basis    9    1   &
basis    10   1   &
basis    11   1   &
basis    11   1   &
basis    13   1   &
basis    14   1   &
basis    15   1   &
basis    16   1

region    bD03  block    INF  INF  ${nydim}  ${nfydim}  INF  INF  units  box
lattice  custom  ${altD}  &
a1  1.0  0.0  0.0  &
a2  0.0  1.0  0.0  &
a3  0.0  0.0  1.0  &
basis  0.0   0.0   0.0   &
basis  0.0   0.5   0.5   &
basis  0.5   0.5   0.0   &
basis  0.5   0.0   0.5   &
basis  0.5   0.0   0.0   &
basis  0.0   0.5   0.0   &
basis  0.0   0.0   0.5   &
basis  0.5   0.5   0.5   &
basis  0.25  0.25  0.25  &
basis  0.75  0.25  0.25  &
basis  0.25  0.75  0.25  &
basis  0.25  0.25  0.75  &
basis  0.25  0.75  0.75  &
basis  0.75  0.25  0.75  &
basis  0.75  0.75  0.25  &
basis  0.75  0.75  0.75  &
orient  x   -1   1   -1  &
orient  y   -1   1    2  &
orient  z    1   1    0
create_atoms  2  region  bD03  &
basis    1    4   &
basis    2    4   &
basis    3    4   &
basis    4    4   &
basis    5    3   &
basis    6    3   &
basis    7    3   &
basis    8    3   &
basis    9    3   &
basis    10   3   &
basis    11   3   &
basis    11   3   &
basis    13   3   &
basis    14   3   &
basis    15   3   &
basis    16   3 

region    rhcp  block   INF INF  -0.05  ${fydim}   INF  INF  units  box
lattice   custom   ${alat} &
a1  1.0    0.0     0.0   &
a2  0.0    ${sqa}  0.0   &
a3  0.0    0.0     ${c}   &
basis   0.0   0.0   0.0   &
basis   0.5   0.5   0.0   &
basis   0.0   ${b1} 0.5   &
basis   0.5   ${b2} 0.5
create_atoms   1  region  rhcp    &
basis   1   5   &
basis   2   5   &
basis   3   5   &
basis   4   5

region    bhcp  block   INF INF  ${fydim} ${ydim}  INF  INF  units  box
lattice   custom   ${alat} &
a1  1.0    0.0     0.0   &
a2  0.0    ${sqa}  0.0   &
a3  0.0    0.0     ${c}   &
basis   0.0   0.0   0.0   &
basis   0.5   0.5   0.0   &
basis   0.0   ${b1} 0.5   &
basis   0.5   ${b2} 0.5
create_atoms   1  region  bhcp    &
basis   1   6   &
basis   2   6   &
basis   3   6   &
basis   4   6

# ---------- Define Interatomic Potential -------------------
pair_style  meam
pair_coeff  *  *   lib_MgNdPb.meam  Mg   Nd  MgNd_para.meam  Mg  Nd  Mg  Nd  Mg  Mg

neighbor      1.0  bin
neigh_modify  every  10  delay  0  check yes

# ---------- Displace atoms and delete overlapping atoms ---------------------
group rhcp region rhcp
group rD03 region rD03
group bhcp region bhcp
group bD03 region bD03
group hcp  union  rhcp  bhcp
group D03  union  rD03  bD03
group  gbregion  union   rhcp  rD03

fix f1  bhcp  setforce 0.0 0.0 0.0
fix f2  bD03  aveforce 0.0 0.0 0.0

delete_atoms     overlap   0.3   rhcp  rD03
displace_atoms   hcp  move  %f  %f  %f  units box

# ---------- Define Settings ---------------------
compute eng all pe/atom
compute eatoms gbregion reduce sum c_eng

thermo   10
thermo_style  custom  step  pe  lx  ly  lz  c_eatoms
dump     mcust   all  custom   200   out/Mg3Nd.*.dump  id  type  mass x y z  c_eng

min_style  cg
minimize   1e-12 1e-10 100000 100000
run        1
variable   gbenergy  equal   "c_eatoms"

print  "gb energy = ${gbenergy}"
                    """ % (param.xdisp,
                           param.ydisp,
                           param.zdisp))
        return
