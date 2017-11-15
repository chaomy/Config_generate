#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-06-28 00:35:14
# @Last Modified by:   chaomy
# @Last Modified time: 2017-11-10 02:08:09


import os
import md_pot_data


class size:
    xdim = None
    ydim = None
    zdim = None


class gb_param:
    xdisp = None
    ydisp = None
    zdisp = None


class gn_md_infile(object):
    def __init__(self, inpot, **kwargs):
        self.pot = inpot
        self.config_file = None
        self.basics = dict([('inrelax', 'in.minimize'),
                            ('innebinit', 'in.init'),
                            ('innebfinal', 'in.final'),
                            ('innpt', 'in.npt'),
                            ('innvt', 'in.nvt'),
                            ('struct', 'bcc'),
                            ('read_data', 'init.txt'),
                            ('read_restart', 'init_restart')])
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                if self.basics.has_key(key):
                    self.basics[key] = value
        return

    def set_md_temperature(self, in_temperature):
        self.md_temperature = in_temperature
        return

    def set_md_config(self, config_file):
        self.config_file = config_file
        return

    def write_cu3au_infile(self,
                           input_size=None):
        from gn_md_input_cu3au import gn_md_input_cu3au
        _drv = gn_md_input_cu3au()

        if input_size is None:
            input_size = size
            input_size.xdim = 70
            input_size.ydim = 50
            input_size.zdim = 20

        _drv._write_Cu3Au_perf("in.test", self.pot["lattice"], input_size)
        return

    def write_mgnd_infile(self, fname, param):
        from gn_md_input_mgnd import gn_md_input_mgnd
        _drv = gn_md_input_mgnd()
        _drv._write_mg_nd_inter_phase(fname, param)
        return

    def write_hcp_lattice_infile(self, fname="in.hcp",
                                 element='Nb'):
        from gn_md_input_hcp_lattice import gn_md_input_hcp_lattice
        _drv = gn_md_input_hcp_lattice()
        _drv._write_hcp_lattice(fname)
        return

    def write_fcc_lattie_infile(self, fname="in.fcc",
                                element='Nb'):
        from gn_md_input_fcc_lattice import gn_md_input_fcc_lattice
        _drv = gn_md_input_fcc_lattice()
        _drv._write_fcc_lattice(fname)
        return

    def write_bcc_lattice_infile(self,
                                 fname="in.bcc"):
        from gn_md_input_bcc_lattice import gn_md_input_bcc_lattice
        _drv = gn_md_input_bcc_lattice()
        _drv._write_bcc_lattice(fname)
        return

    def gn_gsf_minimize(self,
                        config_file="in.gsf", tag='relaxed'):
        from gn_md_input_gsf import gn_md_input_gsf
        _drv = gn_md_input_gsf()
        _drv._write_gsf_minimize(config_file, tag)
        return

    def gn_md_input_vacancy(self,
                            config_file='in.vacancy'):
        from gn_md_input_vacancy import gn_md_input_vacancy
        _drv = gn_md_input_vacancy()
        _drv._write_input_vacancy(config_file)
        return

    def gn_md_input_vacancy_migration(self):
        from gn_md_input_vacancy import gn_md_input_vacancy
        _drv = gn_md_input_vacancy()

        _drv._write_neb_init('init.txt')
        _drv._write_neb_final('final.txt')
        _drv._write_neb_run('init_restart')
        return

    def write_md_thermo_expand(self,
                               *args,
                               **kwargs):
        from gn_md_input_thermo import gn_md_input_thermo
        _drv = gn_md_input_thermo(**self.basics)
        if args is not None:
            for arg in args:
                if arg == 'init':
                    print "find "
                    _drv._write_equilibrium(**kwargs)
                    return
        _drv._write_run_npt(**kwargs)
        return

    def gn_md_minimize_cfg(self,
                           config_file=None, fix=None):

        if config_file is not None:
            self.set_md_config(config_file)
        if fix is None:
            with open("in.minimize", 'w') as fid:
                fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       %s
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom
thermo            1000
thermo_style      custom   step   pe   c_eatoms   etotal
dump         mcfg   all  cfg   100   cfg/%s.*.cfg  mass  type xs ys zs
variable     Ef     equal   "c_eatoms"

#--------------------- MINIMIZE -------------------------------
min_style    cg
minimize     1e-20      1e-20     100000     100000
                        """ % (self.config_file,
                               self.pot["file"],
                               self.pot["element"],
                               self.pot["element"]))
        else:
            with open("in.minimize", 'w') as fid:
                fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       %s
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom
thermo            1000
thermo_style      custom   step   pe   c_eatoms   etotal

group   tofix   id    79  50  51   179  180  181

dump         mcfg   all  cfg   10000   cfg/%s.*.cfg  mass  type xs ys zs
variable     Ef     equal   "c_eatoms"

fix          fixscrew1   tofix    setforce   NULL  NULL  0.0
#--------------------- MINIMIZE -------------------------------
min_style    cg
minimize     1e-20      1e-20     100000     100000
                        """ % (self.config_file,
                               self.pot["file"],
                               self.pot["element"],
                               self.pot["element"]))

        return

    def gn_md_minimize(self,
                       config_file=None):

        if config_file is not None:
            self.set_md_config(config_file)

        with open("in.minimize", 'w') as fid:
            fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       %s
# --------------------- FORCE FIELDS ---------------------
pair_style      %s  
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom
thermo            1000
thermo_style      custom   step   pe   c_eatoms   etotal

variable     Ef     equal   "c_eatoms"

#--------------------- MINIMIZE -------------------------------
min_style    cg
minimize     1e-20      1e-20     100000     100000
run          1
                    """ % (self.config_file,
                           self.pot['pair_style'],
                           self.pot['file'],
                           self.pot['element']))
        return

    def gn_md_shear_lattice(self,
                            config_file=None,
                            temp=None):
        if config_file is not None:
            self.set_md_config(config_file)
        if potential_file is not None:
            self.set_md_potential(potential_file)
        if element is not None:
            self.set_md_element(element)
        with open("in.shear", 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic

# --------------------- FORCE FIELDS ---------------------
read_data       %s
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *   ./%s   %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute     mystress  all   stress/atom     NULL

compute     p1      all     reduce    sum     c_mystress[1]
compute     p2      all     reduce    sum     c_mystress[2]
compute     p3      all     reduce    sum     c_mystress[3]
compute     p4      all     reduce    sum     c_mystress[4]
compute     p5      all     reduce    sum     c_mystress[5]
compute     p6      all     reduce    sum     c_mystress[6]

variable    S11     equal   1e-1*(c_p1)/vol      # Unit is Bar
variable    S22     equal   1e-1*(c_p2)/vol      # we change to Mpa
variable    S33     equal   1e-1*(c_p3)/vol      #
variable    S12     equal   1e-1*(c_p4)/vol      ##xy
variable    S13     equal   1e-1*(c_p5)/vol      ##xz
variable    S23     equal   1e-1*(c_p6)/vol

thermo            10000
thermo_style      custom    step   pe   temp  v_S11  v_S22  v_S33  v_S12  v_S13  v_S23
#--------------------- MINIMIZE -------------------------------
timestep        0.001
variable   srate   equal   1e-4

velocity        all     create     %5.3f   12345   mom     yes    rot     no
fix     1       all     nvt        temp     %5.3f   %5.3f    100
fix     2       all     deform     1    xy  erate   ${srate}   units   box   remap    v

#dump         1       all       cfg    10000   deform_cfg/tensile_*.cfg  mass  type  xs  ys  zs
#dump_modify  1       element   Nb
#dump_modify  1       pad       7
#restart       10000   restart/*.restart

run          500000
unfix   1
unfix   2
                     """ % (self.config_file,
                            self.pot['file'],
                            self.pot["element"],
                            temp, temp, temp))
        return

    def gn_md_temp_lattice(self,
                           config_file=None,
                           temp=None):
        if config_file is not None:
            self.set_md_config(config_file)
        with open("in.lattice", 'w') as fid:
            fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       %s
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom
thermo            1000
thermo_style      custom  step   pe   c_eatoms  lx   ly   lz  pxx  pyy  pzz temp
variable     Ef     equal   "c_eatoms"

#--------------------- MINIMIZE -------------------------------
fix          1    all    box/relax  x  0.0  y  0.0  z  0.0   couple  none  vmax  0.001
min_style    cg
minimize     1e-20      1e-20     100000     100000
unfix        1

timestep        0.001

velocity        all     create     %5.4f   12345   mom     yes    rot     no
fix        2    all     npt        temp     %5.4f   %5.4f    100   &
 iso    0.0   0.0   100  couple  xyz

run       50000

variable lengthx  equal "lx"
variable lengthy  equal "ly"

print "Lattice_constantx = ${lengthx}"
print "Lattice_constanty = ${lengthy}"
                     """ % (self.config_file,
                            self.pot["file"],
                            self.self["element"],
                            temp, temp, temp))
        return

    def gn_md_lattice(self,
                      config_file=None):

        if config_file is not None:
            self.set_md_config(config_file)
        with open("in.lattice", 'w') as fid:
            fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       %s
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *   %s   %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom
thermo            1000
thermo_style      custom  step   pe   c_eatoms   etotal
thermo_style      custom  step   pe   c_eatoms  lx   ly   lz

variable     Ef     equal   "c_eatoms"

#--------------------- MINIMIZE -------------------------------
fix          1    all    box/relax  x  0.0  y  0.0  z  0.0  couple none  vmax 0.001
min_style    cg
minimize     1e-20      1e-20     100000     100000
run          1
unfix        1

variable lengthx  equal "lx"
variable lengthy  equal "ly"

print "Lattice_constantx = ${lengthx}\"
print "Lattice_constanty = ${lengthy}\"
                    """ % (self.config_file,
                           self.pot["file"],
                           self.pot["element"]))
        return

    def gn_md_cij(self):
        with open("init.mod", 'w') as fid:
            fid.write("""
variable up equal 1.0e-6
variable atomjiggle equal 1.0e-5

#units		metal
#variable cfac equal 6.2414e-7
#variable cunits string eV/A^3

units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

variable  etol equal 0.0
variable  ftol equal 1.0e-10
variable  maxiter equal 100
variable  maxeval equal 1000
variable  dmax equal 1.0e-2

variable a equal  %f

boundary	p p p

lattice     bcc $a
region		box prism  0  2.0  0  3.0  0  4.0  0.0  0.0  0.0
create_box	1 box
create_atoms	1 box

mass 1 1.0e-20
                    """ % (self.pot['lattice']))
            fid.close()
        with open("potential.mod", 'w') as fid:
            fid.write("""
pair_style  %s
pair_coeff  * *  %s  %s

neighbor  1.0  nsq
neigh_modify   once  no   every  1  delay  0  check  yes

min_style	     cg
min_modify	     dmax ${dmax} line quadratic

thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
                    """ % (self.pot['pair_style'],
                           self.pot['file'],
                           self.pot['element']))
            fid.close()
        return

    def gn_md_pp_tensile(self,
                         temperature=None,
                         deform_direction=None,
                         deform_rate=None):
        if temperature is not None:
            self.set_md_temperature(temperature)

        with open("in.md_deform", 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       ./relaxed.txt
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *    %s  %s
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
#--------------------- MD Tensile ------------------------------
compute     mystress  all   stress/atom     NULL

compute     p1      all     reduce    sum     c_mystress[1]
compute     p2      all     reduce    sum     c_mystress[2]
compute     p3      all     reduce    sum     c_mystress[3]
compute     p4      all     reduce    sum     c_mystress[4]
compute     p5      all     reduce    sum     c_mystress[5]
compute     p6      all     reduce    sum     c_mystress[6]

variable    S11     equal   1e-1*(c_p1)/vol      # Unit is Bar
variable    S22     equal   1e-1*(c_p2)/vol      # we change to Mpa
variable    S33     equal   1e-1*(c_p3)/vol      #
variable    S12     equal   1e-1*(c_p4)/vol      ##xy
variable    S13     equal   1e-1*(c_p5)/vol      ##xz
variable    S23     equal   1e-1*(c_p6)/vol

reset_timestep	0
timestep        0.001

thermo            1000
thermo_style      custom    step    lx  ly  lz    pe   temp v_S11 v_S22 v_S33 &
v_S12  v_S13  v_S23
#--------------------- MINIMIZE -------------------------------
min_style    cg
minimize     1e-10      1e-10     100000     100000

variable    srate   equal    %s
#--------------------- Temperature -----------------------------
velocity    all     create    %4.3f      12345   mom     yes    rot     no
fix     1   all     nvt       temp    %4.3f      %4.3f      100
fix		2   all     deform    1       %s     erate    ${srate}   units   box   remap    v

variable tmp equal "lx"
variable L0 equal ${tmp}
print "Initial Length, L0: ${L0}"

variable    strain equal   "(lx - v_L0)/ly"
#variable    p1     equal   "v_strain"

fix      def1  all  print   1000   "${strain} ${S11} ${S22} ${S33}"    file  def1.txt  screen  no

dump 		 1      all       cfg   10000   deform_cfg/tensile_*.cfg  mass  type  xs  ys  zs  c_mystress[1]  c_mystress[2] c_mystress[3] c_mystress[4]  c_mystress[5] c_mystress[6]
dump_modify  1      element   %s
dump_modify  1      pad       7

run		 10000000
unfix    1
unfix    2
restart  10000   restart/*.restart
                    """ % (self.potential_file,
                           self.self.pot["element"],
                           deform_rate,
                           2.0 * self.md_temperature,
                           self.md_temperature,
                           self.md_temperature,
                           deform_direction,
                           self.self.pot["element"]))
        return

    def gn_md_nano_tensile(self,
                           temperature=None,
                           deform_direction=None,
                           deform_rate=None):
        if temperature is not None:
            self.set_md_temperature(temperature)

        with open("in.md_deform", 'w') as fid:
            fid.write("""# --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
variable     unit  equal  3.3
# --------------------- FORCE FIELDS ---------------------
read_data       ./lmp_init.txt
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *   %s  %s
neighbor        2.0     bin
neigh_modify    every    1    delay   5  check   yes
# --------------------- SETTINGS ---------------------
compute     peratom     all     pe/atom
compute     eatoms      all     reduce   sum     c_peratom

variable  xsidelow  equal  "xlo"
variable  xsidehi   equal  "xhi"
variable  ysidelo   equal  "ylo"
variable  ysidehi   equal  "yhi"

variable  xrelaxlow  equal (${xsidelow}+3*${unit})
variable  xrelaxhi   equal (${xsidehi}-3*${unit})
variable  yrelaxlow  equal (${ysidelo}+3*${unit})
variable  yrelaxhi   equal (${ysidehi}-3*${unit})

variable  xboundlow  equal (${xsidelow}+${unit})
variable  xboundhi   equal (${xsidehi}-${unit})
variable  yboundlow  equal (${ysidelo}+${unit})
variable  yboundhi   equal (${ysidehi}-${unit})

region   xdeltedlow  block   ${xsidelow}   ${xboundlow}  ${ysidelo}   ${ysidehi}   INF  INF  units box
region   xdeltedhi   block   ${xboundhi}   ${xsidehi}    ${ysidelo}   ${ysidehi}   INF  INF  units  box
region   ydeltedlow  block   ${xsidelow}   ${xsidehi}    ${ysidelo}  ${yboundlow}  INF  INF  units box
region   ydeltedhi   block   ${xsidelow}   ${xsidehi}   ${yboundhi}   ${ysidehi}   INF  INF  units box

region   xfixleft   block   ${xsidelow}   ${xrelaxlow}  ${ysidelo}   ${ysidehi}   INF  INF  units box
region   xfixright  block   ${xrelaxhi}   ${xsidehi}    ${ysidelo}   ${ysidehi}   INF  INF  units  box
region   yfixleft   block   ${xsidelow}   ${xsidehi}    ${ysidelo}  ${yrelaxlow}  INF  INF  units box
region   yfixright  block   ${xsidelow}   ${xsidehi}   ${yrelaxhi}   ${ysidehi}   INF  INF  units box

region   relax    block   ${xrelaxlow}  ${xrelaxhi}  ${yrelaxlow}  ${yrelaxhi}  INF  INF  units  box
#region   fixed    union  4  xfixleft  xfixright  yfixleft  yfixright
region   fixed    union  2  yfixleft  yfixright
region   bound    union  4  xdeltedlow  xdeltedhi  ydeltedlow  ydeltedhi

group    gdelte  region  bound
group    grelax  region  relax
group    gfixed  region  fixed

#fix     a    gfixed      setforce   0.0  0.0  0.0

thermo          10000
thermo_style    custom   step   pe   c_eatoms   etotal
dump            mcfg   all  cfg   10000   cfg/Nb.*.cfg  mass  type xs ys zs

#--------------------- MINIMIZE -------------------------------
min_style    cg
minimize     1e-20      1e-20     100000     100000

compute     mystress  all   stress/atom     NULL

compute     p1      all     reduce    sum     c_mystress[1]
compute     p2      all     reduce    sum     c_mystress[2]
compute     p3      all     reduce    sum     c_mystress[3]
compute     p4      all     reduce    sum     c_mystress[4]
compute     p5      all     reduce    sum     c_mystress[5]
compute     p6      all     reduce    sum     c_mystress[6]

variable    S11     equal   1e-1*(c_p1)/vol      # Unit is Bar
variable    S22     equal   1e-1*(c_p2)/vol      # we change to Mpa
variable    S33     equal   1e-1*(c_p3)/vol      #
variable    S12     equal   1e-1*(c_p4)/vol      ##xy
variable    S13     equal   1e-1*(c_p5)/vol      ##xz
variable    S23     equal   1e-1*(c_p6)/vol

undump          mcfg
timestep        0.001

thermo            10000
thermo_style      custom    step    lx  ly  lz   pe   temp

variable    srate   equal    %s
#--------------------- Temperature -----------------------------
velocity    all     create    %s      12345   mom     yes    rot     no
fix     1   all     nvt       temp    %s      %s      100
fix		2   all     deform    1       %s      erate    ${srate}   units   box   remap    v

variable tmp equal "lx"
variable L0 equal ${tmp}
print "Initial Length, L0: ${L0}"

variable    strain equal   "(lx - v_L0)/ly"
#variable    p1     equal   "v_strain"
#variable    p2     equal   "-pxx/10000"
#variable    p3     equal   "-pyy/10000"
#variable    p4     equal   "-pzz/10000"

fix      def1  all  print   10000   "${strain} ${S11} ${S22} ${S33} ${S12} ${S13} ${S23}"   file  def1.txt  screen  no

dump 		 1      all       cfg   10000   cfg/tensile_*.cfg  mass  type  xs  ys  zs
dump_modify  1      element   %s
dump_modify  1      pad       7

run		5000000
unfix   1
unfix   2

restart  10000  restart/*.restart
                    """ % (self.potential_file,
                           self.self.pot["element"],
                           deform_rate,
                           self.md_temperature,
                           self.md_temperature,
                           self.md_temperature,
                           deform_direction,
                           self.self.pot["element"]))
            fid.close()
            return

    def gn_md_add_force(self,
                        temp,
                        stress):
        with open("in.md_addforce", 'w') as fid:
            fid.write("""
                    # --------------------- INITIALIZAITION ---------------------
clear
units	     metal
dimension	 3
boundary	 p  p  p
atom_style	 atomic
# --------------------- FORCE FIELDS ---------------------
read_data       ./relaxed.txt
# --------------------- FORCE FIELDS ---------------------
pair_style      eam/alloy
pair_coeff      *  *    W.set.txt   W
neighbor        2.0     bin
neigh_modify    every    1    delay   0  check   yes
#--------------------- MD Tensile ------------------------------
compute     mystress  all   stress/atom     NULL

variable  latticeC  equal  3.1648492
variable  ysidelo   equal  "ylo"
variable  ysidehi   equal  "yhi"

compute     p1      all     reduce    sum     c_mystress[1]
compute     p2      all     reduce    sum     c_mystress[2]
compute     p3      all     reduce    sum     c_mystress[3]
compute     p4      all     reduce    sum     c_mystress[4]
compute     p5      all     reduce    sum     c_mystress[5]
compute     p6      all     reduce    sum     c_mystress[6]

variable    S11     equal   1e-1*(c_p1)/(0.75*vol)      # Unit is Bar
variable    S22     equal   1e-1*(c_p2)/(0.75*vol)      # we change to Mpa
variable    S33     equal   1e-1*(c_p3)/(0.75*vol)      #
variable    S12     equal   1e-1*(c_p4)/(0.75*vol)      ##xy
variable    S13     equal   1e-1*(c_p5)/(0.75*vol)      ##xz
variable    S23     equal   1e-1*(c_p6)/(0.75*vol)

#variable   ytop    equal  (${ysidehi}-(11.332813*${latticeC}))
#variable   ysurfhi  equal   (${ysidehi}-(11.332813*${latticeC}))
#variable   ysurflo  equal   (${ysidelo}+(11.0265240*${latticeC}))

variable   ysurfhi  equal   (${ysidehi}-41.34)
variable   ysurflo  equal   (${ysidelo}+41.34)

#region   surftop   block    INF  INF   ${ytop}     ${ysidehi}  INF  INF units box
region   surfup    block    INF  INF    ${ysurfhi}  ${ysidehi}  INF  INF units box
region   surfdown  block    INF  INF    ${ysidelo}  ${ysurflo}  INF  INF units box
region   relax     block    INF  INF    ${ysurflo}  ${ysurfhi}  INF  INF units box

#group    gtop       region  surftop
#delete_atoms    group   gtop

group    gsurfup    region  surfup
group    gsurfdown  region  surfdown
group    relax      region  relax

reset_timestep	0
timestep        0.001

thermo            2000
thermo_style      custom    step    lx  ly  lz    pe   temp  v_S11  v_S22  v_S33  v_S12  v_S13  v_S23

#--------------------- Temperature -----------------------------
variable    mytemp   equal   %d
variable    inittemp equal   2*${mytemp}

velocity    relax      create   ${inittemp}      12345   mom     yes   rot     no

velocity    gsurfup    create   000.000      12345   mom     yes    rot     no
velocity    gsurfdown  create   000.000      12345   mom     yes    rot     no

variable    shearforceup   equal        %s
variable    shearforcedown equal       -%s

fix     1   all     nvt     temp    ${mytemp}    ${mytemp}    100
fix     2   gsurfup     addforce    0.0   0.0   ${shearforceup}
fix     3   gsurfdown   addforce    0.0   0.0   ${shearforcedown}

dump 		 1      all       cfg   10000   force_cfg/tensile_*.cfg  mass  type  xs  ys  zs
dump_modify  1      element   Nb
dump_modify  1      pad       7
restart      100000   restart/*.restart

run		 10000000

unfix    1
unfix    2
unfix    3
                    """ % (temp,
                           stress,
                           stress))
        return

    def gn_md_tensile(self):
        with open("in.stat_tensile", 'w') as fid:
            fid.write("""
# --------------------- INITIALIZAITION ---------------------
clear
units		metal
dimension	3
boundary	p p p
atom_style	atomic

# --------------------- ATOM DEFINITION ---------------------
read_data   ./lattice.txt

variable    xlen    equal   lx
variable    ylen    equal   ly
variable    zlen    equal   lz

# --------------------- FORCE FIELDS ---------------------
pair_style  eam/alloy
pair_coeff * *   %s  %s
neighbor        2.0     bin
neigh_modify    every   1     delay   0   check   yes
# --------------------- SETTINGS ---------------------
timestep     0.001     ##### 10e-15

compute      mystress  all    stress/atom   NULL

compute     p1      all     reduce sum c_mystress[1]
compute     p2      all     reduce sum c_mystress[2]
compute     p3      all     reduce sum c_mystress[3]
compute     p4      all     reduce sum c_mystress[4]
compute     p5      all     reduce sum c_mystress[5]
compute     p6      all     reduce sum c_mystress[6]

variable         S11 equal (c_p1)/vol       ##xx
variable         S22 equal (c_p2)/vol
variable         S33 equal (c_p3)/vol
variable         S12 equal (c_p4)/vol       ##xy
variable         S13 equal (c_p5)/vol       ##xz
variable         S23 equal (c_p6)/vol       ##yz

#----------------------------------------------------
thermo   100
thermo_style  custom  step  pe   lx   ly   lz   v_S11  v_S22  v_S33  v_S12  v_S13  v_S23
#--------------------relax box -----------------------
min_style   cg
minimize    10e-20 10e-20  200000   200000

variable Sxx equal "v_S11/10000"
variable Syy equal "v_S22/10000"
variable Szz equal "v_S33/10000"
variable Sxy equal "v_S12/10000"
variable Sxz equal "v_S13/10000"
variable Syz equal "v_S23/10000"

####################################
print "Lx = ${xlen}"
print "Ly = ${ylen}"
print "Lz = ${zlen}"

print "Sxx = ${Sxx}"
print "Syy = ${Syy}"
print "Szz = ${Szz}"
print "Syz = ${Syz}"
print "Sxz = ${Sxz}"
                    """ % (self.potential_file, self.self.pot["element"]))
            fid.close()
            return

    def write_shear_infile(self, file, basis, box):
        box = box * 10
        with open("in.tmp", 'w') as fid:
            fid.write("""
                # ---------- Initialize Simulation ---------------------
clear
units       metal
atom_style  atomic
boundary    p   p   p
atom_modify sort  0 0.0

# ---------- Create Atomistic Structure ------------
variable  alat  equal  3.308

lattice  custom    ${alat}   &
a1       %f   %f   %f    &
a2       %f   %f   %f    &
a3       %f   %f   %f    &
basis    0.0   0.0   0.0  &
basis    0.5   0.5   0.5

region     whole   prism     0    %f    0     %f    0     %f    %f    %f    %f    units box
create_box     2       whole
create_atoms   2       region    whole  &
basis    1     1      &
basis    2     2

pair_style  eam/alloy
pair_coeff  *  *  ../Nb.eam.alloy.webarchive   Nb  Nb
neighbor   1.0  bin
neigh_modify  every  1  delay 0 check yes

dump       mcust   all   custom   500    tmp/%s  id  type  mass x y z

minimize   1e-8   1e-8   100000   100000
                    """ % (
                basis[0, 0], basis[0, 1], basis[0, 2],
                basis[1, 0], basis[1, 1], basis[1, 2],
                basis[2, 0], basis[2, 1], basis[2, 2],
                box[0, 0],  box[1, 1], box[2, 2],
                box[1, 0], box[2, 0], box[2, 1],
                file))
        return


if __name__ == '__main__':
    drv = gn_md_infile()
    drv.gn_md_cij()