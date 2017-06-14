#!/usr/local/bin/python
# encoding: utf-8
#
import glob
import numpy as np
import matplotlib.pylab as plt
from   mpl_toolkits.mplot3d import Axes3D
from   mpi4py import MPI
import  sys, os

class MD_ChangeBox(object):
    def __init__(self):
        self.lattice_constant = 3.307898
        self.pi =  3.141592653589793
        self.burger = np.sqrt(3)/2. * self.lattice_constant
        self.screw_coeff = self.burger /(2. * self.pi)
        self.Edge_coeff = self.burger/(2. * self.pi)
        self.P = 0.33
        return

    def read_xyz(self,filename):
        xs,ys,zs,Atomid  = [],[],[],[]
        with open(filename, 'r') as fid :
            Raw = fid.readlines()
        fid.close()
        AtomNumber = int(Raw[3])
        xlo , xhi = float(Raw[5].split()[0]) , float(Raw[5].split()[1])
        ylo , yhi = float(Raw[6].split()[0]) , float(Raw[6].split()[1])
        zlo , zhi = float(Raw[7].split()[0]) , float(Raw[7].split()[1])
        for i in range(9,AtomNumber + 9):
            xs.append(float(Raw[i].split()[1]))
            ys.append(float(Raw[i].split()[2]))
            zs.append(float(Raw[i].split()[3]))
            Atomid.append(float(Raw[i].split()[-1]))
        return np.array(xs), np.array(ys), np.array(zs) , AtomNumber, xlo, xhi, ylo, yhi, zlo, zhi, np.array(Atomid)

    def read_cfg(self,filename):
        import re
        x,y,z,energy,ccsym,Sxx,Syy,Szz,Sxy,Sxz,Syz,coord1,atomid = [],[],[],[],[],[],[],[],[],[],[],[],[]
        with open(filename, 'r') as fid:
            for line in fid :
                match = re.search(r"^([+\-]?(\d+(\.\d*)?|\d*\.\d+)([eE][+\-]?\d+)?) ([+\-]?(\d+(\.\d*)?|\d*\.\d+)([eE][+\-]?\d+)?)" ,line)
                if match: # xs ys zs  c_csym c_peratom c_mystress[1] c_mystress[2] c_mystress[3] c_mystress[4] c_mystress[5] c_mystress[6] id  c_coordn
                    raw = line.split()
                    x.append(float(raw[0]))
                    y.append(float(raw[1]))
                    z.append(float(raw[2]))
                    ccsym.append(float(raw[3]))
                    energy.append(float(raw[4]))
                    Sxx.append(float(raw[5]))
                    Syy.append(float(raw[6]))
                    Szz.append(float(raw[7]))
                    Sxy.append(float(raw[8]))
                    Sxz.append(float(raw[9]))
                    Syz.append(float(raw[10]))
                    atomid.append(int(raw[11]))
                    coord1.append(int(raw[12]))
                match2 = re.search(r"H0\(1,1\)\s*=\s*(\d*.\d*)",line)
                if match2 :
                    Lx = float(match2.group(1))
                match3 = re.search(r"H0\(2,2\)\s*=\s*(\d*.\d*)",line)
                if match3 :
                    Ly = float(match3.group(1))
                match4 = re.search(r"H0\(3,3\)\s*=\s*(\d*.\d*)",line)
                if match4:
                    Lz = float(match4.group(1))
            fid.close()
        return x,y,z,ccsym,energy,Sxx,Syy,Szz,Sxy,Sxz,Syz,coord1,atomid,Lx,Ly,Lz

    def Intro_Screw_cfg(self):
        filename = sys.argv[1]
        lattice_constant = 3.30789893315
        x,y,z,ccsym,energy,Sxx,Syy,Szz,Sxy,Sxz,Syz,coord1,atomid,Lx0,Ly0,Lz0 = self.read_cfg(filename)
        xnew,ynew,znew,coordN,energyN,atomidN = [],[],[],[],[],[]
        a =   float(Ly0)/ float(lattice_constant)
        dd = 0.4 /a
        xc ,yc = 0.5,0.5
        AtomNumber = len(x)
        for i in range(AtomNumber):
            if x[i] > xc :
                if  np.abs(y[i] - yc) > dd :
                    xnew.append(float(x[i])) ; ynew.append(float(y[i])); znew.append(float(z[i]))
                    coordN.append(coord1[i]) ; atomidN.append(atomid[i]); energyN.append(energy[i])
            else :
                xnew.append(float(x[i])) ; ynew.append(float(y[i])); znew.append(float(z[i]))
                coordN.append(coord1[i]) ; atomidN.append(atomid[i]); energyN.append(energy[i])
        print "atom total: " , AtomNumber
        print "atom deleted" , AtomNumber - len(xnew)
        coeff = self.screw_coeff
        for i in range(len(xnew)):
            theta = np.arctan2((ynew[i] - yc),(xnew[i] - xc))
            dz = self.lattice_constant * coeff * theta
            print dz
            znew[i] += dz
        with open("screw.cfg",'w') as fid:
            fid.write("Number of particles = %d\n"%(len(xnew)))
            fid.write("A = 1 Angstrom (basic length-scale)\n")
            fid.write("""H0(1,1) = %f A
H0(1,2) = 0 A
H0(1,3) = 0 A
H0(2,1) = 0 A
H0(2,2) = %f A
H0(2,3) = 0 A
H0(3,1) = 0 A
H0(3,2) = 0 A
H0(3,3) = %f A
.NO_VELOCITY.
entry_count = 6
auxiliary[0] = c_coordn
auxiliary[1] = id
auxiliary[2] = energy
"""%(Lx0,Ly0,Lz0))
            for i in range(len(xnew)):
                fid.write("""92.906380
C
%6.4f %6.4f %6.4f %d %d %6.4f
"""%(xnew[i],ynew[i],znew[i],coordN[i],atomidN[i],energyN[i]))
        return

    def Intro_Screw_dipole(self):
        filelist = glob.glob('./xyz/*')
        filename = filelist[-1]
        xs , ys, zs , AtomNumber ,xlo, xhi, ylo, yhi, zlo, zhi, Atomid = self.read_xyz(filename)
        xs = np.array(xs) ; ys = np.array(ys) ; zs = np.array(zs)
        xc1 , xc2 , yc = 0.25 * xhi + 1.25 * xlo , 0.75 * xhi - 0.25 * xlo , 0.5 * (yhi + ylo)
        #-------------------- cut  plane along x positive direction
        coeff  = self.screw_coeff
        xnew , ynew , znew = [] , [], []
        for i in range(AtomNumber):
            xnew.append(float(xs[i])) ; ynew.append(float(ys[i])); znew.append(float(zs[i]))
        AtomNumberOut = len(xnew)
        print AtomNumber
        print "delete %d atoms"%(AtomNumber-AtomNumberOut)
        for i in range(AtomNumberOut):
            theta = np.arctan2((ynew[i] - yc),(xnew[i] - xc1))
            dz = coeff * theta
            znew[i] = dz + znew[i]
        for i in range(AtomNumberOut):
            theta = np.arctan2((ynew[i] - yc),(xnew[i] - xc2))
            dz = -coeff * theta
            znew[i] = dz + znew[i]
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_zlim(-2.5,6.5)
        ax.scatter(xnew,ynew,znew)
        plt.show()

        filename = "screw.txt"
        with open(filename,mode="w") as fout:
            fout.write("#screw_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(Atomid[i], xnew[i] ,ynew[i] ,znew[i] ))
        fout.close()
        return

    def Mpi_Intro_Screw(self):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        def update_position(x, y, z, xc, yc, coeff):
            for i in range(len(x)):
                dx, dy = x[i] - xc , y[i] - yc
                theta = np.arctan2(dy, dx)
                dz = coeff * theta
                z[i] += dz
            return (x, y, z)
        if rank == 0 :
            filelist = glob.glob("./xyz/*")
            filename = filelist[-1]
            xs , ys, zs , AtomNumber ,xlo, xhi, ylo, yhi, zlo, zhi, Atomid = self.read_xyz(filename)
            xc = 0.5 * (xlo + xhi) ; yc = 0.5 * (ylo + yhi)
            local_n = AtomNumber/size
            LocalN_xc_yc = np.concatenate([local_n, xc, yc])
        else :
            LocalN_xc_yc  = np.array([0])
            x, y, z = None, None, None
        comm.Bcast(LocalN_xc_yc, root = 0)
        coeff  = self.screw_coeff

        if rank == 0 :
            if (len(x) % size == 0):
                check = True
            else :
                check = False
                body = (len(x) - np.mod(len(x), size))
                xrem, yrem, zrem = x[body:], y[body:], z[body:]
                print "Now in rank 0  x  and y size" , len(x), len(y) , len(z)
        x_local, y_local, z_local  = np.zeros(LocalN_xc_yc[0]), np.zeros(LocalN_xc_yc[0]) , np.zeros(LocalN_xc_yc[0])

        comm.Scatter(x, x_local, root=0)
        comm.Scatter(y, y_local, root=0)
        comm.Scatter(z, z_local, root=0)

        (x_local, y_local, z_local) = update_position(x_local, y_local, z_local, LocalN_xc_yc[1], LocalN_xc_yc[2], coeff)

        if rank == 0 :
            if check == False :
                (xrem , yrem, zrem) = update_position(xrem, yrem, zrem, LocalN_xc_yc[1], LocalN_xc_yc[2], coeff)

        comm.Gather(z_local, z, root = 0)
        if rank == 0 :
            if check == False :
                x, y, z = np.concatenate([x, xrem]), np.concatenate([y, yrem]), np.concatenate([z, zrem])
            filename = "Screw.txt"
            with open(filename,mode="w") as fout:
                fout.write("#Edge_dislocation for bcc [111](110)\n")
                fout.write("\n")
                fout.write("%d atoms\n"%(AtomNumber))
                fout.write("1 atom types\n")
                fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
                fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
                fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
                fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
                fout.write("\n")
                fout.write("Atoms\n")
                fout.write("\n")
                for i in range(AtomNumber):
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(Atomid[i], x[i],y[i] ,z[i]))
                fout.close()

        return

    def Intro_Screw(self):
        filelist = glob.glob("./xyz/*")
        filename = filelist[-1]

        xs , ys, zs , AtomNumber ,xlo, xhi, ylo, yhi, zlo, zhi, Atomid = self.read_xyz(filename)
        xs = np.array(xs) ; ys = np.array(ys) ; zs = np.array(zs)
        xc = 0.5 * (xlo + xhi) ; yc = 0.5 * (ylo + yhi)
        #-------------------- cut  plane along x positive direction
        coeff  = self.screw_coeff
        xnew , ynew , znew = [] , [], []
        for i in range(AtomNumber):
            xnew.append(float(xs[i])) ; ynew.append(float(ys[i])); znew.append(float(zs[i]))
            # ------------------- intro screw dislocation -------------
        AtomNumberOut = len(xnew)
        print AtomNumber
        for i in range(AtomNumberOut):
            dx, dy = xnew[i] - xc , ynew[i] - yc
            theta = np.arctan2(dy, dx)
            dz = coeff * theta
            znew[i] = dz + znew[i]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        #ax.set_zlim(-2.5,10)
        ax.scatter(xnew,ynew,znew)
        plt.show()
        # ------------------- OutPut screw.txt --------------------
        filename = "screw.txt"
        with open(filename,mode="w") as fout:
            fout.write("#screw_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(Atomid, xnew[i] ,ynew[i] ,znew[i] ))
        fout.close()
        return

    def Intro_Edge(self):
        import math
        filelist = glob.glob("xyz/*")
        filename = filelist[-1]
        xs , ys, zs , AtomNumber ,xlo, xhi, ylo, yhi, zlo, zhi, Atomid  = self.read_xyz(filename)

        xs = np.array(xs) ; ys = np.array(ys) ; zs = np.array(zs)
        xc = 0.5 * (np.max(xs) + np.min(xs))
        yc = 0.5 * (np.max(ys) + np.min(ys))
        #-------------------- cut  plane along x positive direction
        coeff  = self.Edge_coeff
            # ----------------- intro Edge dislocation -------------
        AtomNumberOut = len(xs)
        print   "Atom Number ", AtomNumber

        for i in range(AtomNumberOut):
            dx , dy = xs[i] - xc , ys[i] - yc
            A = (1 - self.P) * (dx**2 + dy**2)
            ux = coeff * (math.atan2(dy , dx) + (dx * dy)/(2 * A))
            uy = (1 -2*self.P)/(4 *(1 - self.P)) * math.log(dx**2 + dy**2) + (dx**2-dy**2)/(4 * (1 - self.P)*(dx**2 + dy**2))
            uy *= -coeff
            xs[i] += 0.1 *  ux
            ys[i] += 0.1 *  uy

        xlo , xhi , ylo, yhi, zlo, zhi = min(xs) - 4  , max(xs) + 4, min(ys) - 4, max(yhi) + 4, min(zs) -4 , max(zs) + 4

        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.set_zlim(-2.5,6.5)
        ax.scatter(xs,ys,zs)
        plt.show()
        # ------------------- OutPut Edge.txt --------------------
        filename = "Edge.txt"
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(Atomid[i], xs[i] ,ys[i] ,zs[i] ))
        fout.close()
        return

    def Intro_Crack_cfg(self):
        import math
        filename = sys.argv[1]
        lattice_constant = 3.30789893315
        x,y,z,ccsym,energy,Sxx,Syy,Szz,Sxy,Sxz,Syz,coord1,atomid,Lx0,Ly0,Lz0 = self.read_cfg(filename)
        xnew,ynew,znew,coordN,energyN,atomidN = [],[],[],[],[],[]
        a = lattice_constant / Ly0
        pi = self.pi
        xc ,yc = 0.2 ,0.5
        AtomNumber = len(x)
        for i in range(AtomNumber):
            if z[i] < 0.3 and z[i] > -0.01:
                if x[i] < xc and y[i] - yc < 0.01:
                    xnew.append(float(x[i])) ; ynew.append(float(y[i])); znew.append(float(z[i]))
                    coordN.append(coord1[i]) ; atomidN.append(atomid[i]); energyN.append(energy[i])
                    print "delete"
                else :
                    xnew.append(float(x[i])) ; ynew.append(float(y[i])); znew.append(float(z[i]))
                    coordN.append(coord1[i]) ; atomidN.append(atomid[i]); energyN.append(energy[i])
        stress = 5.0 ; k1 = stress * np.sqrt(self.pi * lattice_constant) ; Pratio = 0.333 ; shearModule = 50.0
        coeff = k1 * np.sqrt(self.pi * 2)/(8 * shearModule * self.pi)
        print coeff
        for i in range(len(xnew)):
            dx , dy = xnew[i] - xc , ynew[i] - yc
            if dx == 0 : dx = 0.000000000001
            if dy == 0 : dy = 0.000000000001
            r = math.sqrt(dx*dx + dy*dy)
            if ynew[i] > yc:
                if xnew[i] > xc :
                    theta = math.atan((dy/dx))
                if xnew[i] < xc :
                    theta = math.atan((-dx/dy)) + 0.5 * pi
                ux = coeff * math.sqrt(r) * ((5-8*Pratio)*math.cos(0.5 * theta) - math.cos(1.5 * theta))
                uy = coeff * math.sqrt(r) * ((7-8*Pratio)*math.sin(0.5 * theta) - math.sin(1.5 * theta))
            if ynew[i] < yc:
                if xnew[i] > xc :
                    theta = math.atan((-dy/dx))
                if xnew[i] < xc :
                    theta = math.atan((dx/dy)) + 0.5 * pi
                ux = coeff * math.sqrt(r) * ((5-8*Pratio)*math.cos(0.5 * theta) - math.cos(1.5 * theta))
                uy = -1 *  coeff * math.sqrt(r) * ((7-8*Pratio)*math.sin(0.5 * theta) - math.sin(1.5 * theta))
            xnew[i] += ux
            ynew[i] += uy
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection='3d')
    #    ax.set_zlim(-2.5,6.5)
        ax.scatter(xnew,ynew,znew)
        plt.show()
        AtomNumberOut = len(xnew)
    	filename = "Crack.txt"
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumberOut))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(0,Lx0))
            fout.write("%f\t%f ylo yhi\n"%(0,Ly0))
            fout.write("%f\t%f zlo zhi\n"%(0,Lz0))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumberOut):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(i+1, xnew[i],ynew[i] ,znew[i] ))
        fout.close()
        return

    def Mpi_Intro_Crack_xyz(self,k1,p1,p2,q1,q2,u1,u2):
        from numpy.lib import scimath as SM
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        def updata_position(x, y, xc, yc):
            for i in range(len(x)):
                dx , dy = x[i] - xc , y[i] - yc
                r = np.sqrt(dx * dx + dy * dy)
                coeff = 100 * k1 * np.sqrt(2 * r / self.pi)
                if y[i] < yc:
                    theta = np.arctan2(-dy, dx)
                    ux  = 1./(u1-u2) * (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    uy  = 1./(u1-u2) * (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux  = coeff * ux.real
                    uy  = coeff * uy.real
                if y[i] >= yc:
                    theta = np.arctan2(dy, dx)
                    ux  = 1./(u1-u2) * (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    uy  = 1./(u1-u2) * (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                    ux  = coeff * ux.real
                    uy  = -coeff * uy.real
                print ux
                print uy
                x[i] += ux
                y[i] += uy
            return (x, y)

        if rank == 0 :
            filelist = glob.glob("./xyz/*")
            filename = filelist[-1]
            x, y, z , AtomNumber, xlo, xhi, ylo, yhi, zlo, zhi, Atomid  = self.read_xyz(filename)
            local_n = np.array([AtomNumber/size])
            print "in rank 0  x  and y size" , len(x), len(y)
        else :
            local_n = np.array([0])
            x, y = None, None

        comm.Bcast(local_n, root = 0)
        print "local_n of %d is "%(rank), local_n

        if rank == 0:
            if (len(x) % size == 0):
                check = True
            else :
                check = False
                body = (len(x) - np.mod(len(x),size))
                xrem, yrem = x[body:], y[body:]
                x , y = x[0:body], y[0:body]
                xrem , yrem = np.array(xrem), np.array(yrem)
            print "Now in rank 0  x  and y size" , len(x), len(y)
        else :
            print "I am rank %d"%(rank)

        xc , yc = 0, 0
        x_local = np.zeros(local_n)
        y_local = np.zeros(local_n)

        if rank == 0 :
            print rank,  len(x_local) , len(y_local), len(x), len(y)
        if rank != 0 :
            print rank,  len(x_local) , len(y_local), x, y

        comm.Scatter(x, x_local, root=0)
        comm.Scatter(y, y_local, root=0)

        (x_local, y_local)  = updata_position(x_local, y_local, xc, yc)

        if rank == 0 :
            if check == False :
                (xrem, yrem) = updata_position(xrem, yrem, xc, yc)

        comm.Gather(x_local, x, root = 0)
        comm.Gather(y_local, y, root = 0)

        if rank == 0 :
            if check == False :
                x = np.concatenate([x, xrem])
                y = np.concatenate([y, yrem])
            filename = "Crack.txt"
            xlo , xhi, ylo, yhi = min(x) - 4 , max(x) + 4 , min(y) - 4 , max(y) + 4
            with open(filename,mode="w") as fout:
                fout.write("#Edge_dislocation for bcc [111](110)\n")
                fout.write("\n")
                fout.write("%d atoms\n"%(AtomNumber))
                fout.write("1 atom types\n")
                fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
                fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
                fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
                fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
                fout.write("\n")
                fout.write("Atoms\n")
                fout.write("\n")
                for i in range(AtomNumber):
                    fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(Atomid[i], x[i],y[i] ,z[i]))
                fout.close()
        return

    def Intro_Crack_xyz(self,k1,p1,p2,q1,q2,u1,u2):
        from numpy.lib import scimath as SM
        filelist = glob.glob("./xyz/*")
        filename = filelist[-1]
        x, y, z , AtomNumber, xlo, xhi, ylo, yhi, zlo, zhi, Atomid  = self.read_xyz(filename)
        xc ,yc = 0 ,0
        AtomNumber = len(x)
        for i in range(len(x)):
            dx , dy = x[i] - xc , y[i] - yc
            r = np.sqrt(dx * dx + dy * dy)
            coeff = 100 * k1 * np.sqrt(2 * r / self.pi)
            if y[i] < yc:
                theta = np.arctan2(-dy, dx)
                ux  = 1./(u1-u2) * (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                uy  = 1./(u1-u2) * (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                ux  = coeff * ux.real
                uy  = coeff * uy.real
            if y[i] >= yc:
                theta = np.arctan2(dy, dx)
                ux  = 1./(u1-u2) * (u1 * p2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * p1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                uy  = 1./(u1-u2) * (u1 * q2 * SM.sqrt(np.cos(theta) - u2 * np.sin(theta)) - u2 * q1 * SM.sqrt(np.cos(theta) - u1 * np.sin(theta)))
                ux  = coeff * ux.real
                uy  = -coeff * uy.real
            print ux
            print uy
            x[i] += ux
            y[i] += uy
        xlo , xhi, ylo, yhi = min(x) - 4 , max(x) + 4 , min(y) - 4 , max(y) + 4
   #     import matplotlib.pylab as plt
   #     from mpl_toolkits.mplot3d import Axes3D
   #     fig = plt.figure()
   #     ax = fig.gca(projection='3d')
   #     ax.scatter(xnew,ynew,znew)
   #     plt.show()
    	filename = "Crack.txt"
        with open(filename,mode="w") as fout:
            fout.write("#Edge_dislocation for bcc [111](110)\n")
            fout.write("\n")
            fout.write("%d atoms\n"%(AtomNumber))
            fout.write("1 atom types\n")
            fout.write("%f\t%f xlo xhi\n"%(xlo,xhi))
            fout.write("%f\t%f ylo yhi\n"%(ylo,yhi))
            fout.write("%f\t%f zlo zhi\n"%(zlo,zhi))
            fout.write("xy xz yz = %8.5f %8.5f %8.5f"%(0,0,0))
            fout.write("\n")
            fout.write("Atoms\n")
            fout.write("\n")
            for i in range(AtomNumber):
                fout.write("%d  1  %12.7f %12.7f %12.7f\n"%(Atomid[i], xnew[i],ynew[i] ,znew[i]))
        fout.close()
        return

class VA_ChangeBox(object):
    def __init__(self):
        self.lattice_constant = 3.3078
        self.pi = 3.141592653589793
        self.burger = np.sqrt(3)/2. * self.lattice_constant
        self.screw_coeff = self.burger/(2. * self.pi)
        self.Edge_coeff = self.burger /(2. * self.pi)
        return
    def read_pos(self):
        with open("POSCAR",'r')  as fid :
            Raw = fid.readlines()
            fid.close()
        self.lattice_constant = float(Raw[1])
        print self.lattice_constant
        H = np.zeros([3,3],"float")
        for i in range(3):
            for j in range(3):
                H[i,j] = float(Raw[2 + i].split()[j])
        AtomNumber = int(Raw[5])
        print AtomNumber
        Comm = str(Raw[6])
        print Comm
        AtomPosition = np.zeros([3, AtomNumber],"float")
        for j in range(AtomNumber):
            for i in range(3):
                AtomPosition[i,j] = float(Raw[7 + j].split()[i])
        print H
        print AtomPosition
        return AtomNumber, np.mat(H), Comm, AtomPosition

    def volume_conserving_Ortho_strain(self,delta):
        AtomNumber, H, Comm, AtomPosition = self.read_pos()
        OrthoM = np.mat([[1 + delta, 0,0],[0,1 - delta,0],[0,0,1/(1 - delta**2)]])
        Transformed_Base = OrthoM * H
        with open("POSCAR", 'w') as fid:
            fid.write("# Screw Bcc\n")
            fid.write("%12.6f\n"%(self.lattice_constant))
            fid.write("%12.6f %12.6f %12.6f\n"%(Transformed_Base[0,0],Transformed_Base[0,1],Transformed_Base[0,2]))
            fid.write("%12.6f %12.6f %12.6f\n"%(Transformed_Base[1,0],Transformed_Base[1,1],Transformed_Base[1,2]))
            fid.write("%12.6f %12.6f %12.6f\n"%(Transformed_Base[2,0],Transformed_Base[2,1],Transformed_Base[2,2]))
            fid.write("%d\n"%(AtomNumber))
            fid.write(Comm)
            for i in range(AtomNumber):
                fid.write("%12.6f %12.6f %12.6f\n"%(AtomPosition[0,i],AtomPosition[1,i],AtomPosition[2,i]))
            fid.close()
        return

    def volume_conserving_Mono_strain(self, delta):
        AtomNumber, H, Comm, AtomPosition = self.read_pos()
        OrthoM = np.mat([[1 ,delta, 0],[delta ,1 ,0],[0,0,1/(1 - delta**2)]])
        Transformed_Base = OrthoM * H
        with open("POSCAR", 'w') as fid:
            fid.write("# Screw Bcc\n")
            fid.write("%12.6f\n"%(self.lattice_constant))
            fid.write("%12.6f %12.6f %12.6f\n"%(Transformed_Base[0,0],Transformed_Base[0,1],Transformed_Base[0,2]))
            fid.write("%12.6f %12.6f %12.6f\n"%(Transformed_Base[1,0],Transformed_Base[1,1],Transformed_Base[1,2]))
            fid.write("%12.6f %12.6f %12.6f\n"%(Transformed_Base[2,0],Transformed_Base[2,1],Transformed_Base[2,2]))
            fid.write("%d\n"%(AtomNumber))
            fid.write(Comm)
            for i in range(AtomNumber):
                fid.write("%12.6f %12.6f %12.6f\n"%(AtomPosition[0,i],AtomPosition[1,i],AtomPosition[2,i]))
            fid.close()
        return

    def Intro_Screw(self):
        AtomNumber, H, Comm, AtomPosition = self.read_pos()
        xc = 0.5 * H[0,0] ; yc = 0.5 * H[1,1]
        print xc , yc
        for i in range(AtomNumber):
            dx , dy = AtomPosition[0, i] - xc , AtomPosition[1, i] - yc
            theta = np.arctan2(dy, dx)
            dz = self.screw_coeff * theta
            AtomPosition[2, i] += dz
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.scatter(AtomPosition[0,:],AtomPosition[1,:],AtomPosition[2,:])
        plt.show()
        with open("POSCAR", 'w')  as fid:
            fid.write("# Screw Bcc\n")
            fid.write("%12.6f\n"%(self.lattice_constant))
            fid.write("%12.6f %12.6f %12.6f\n"%(H[0,0],H[0,1],H[0,2]))
            fid.write("%12.6f %12.6f %12.6f\n"%(H[1,0],H[1,1],H[1,2]))
            fid.write("%12.6f %12.6f %12.6f\n"%(H[2,0],H[2,1],H[2,2]))
            fid.write("%d\n"%(AtomNumber))
            fid.write("Cartesian\n")
            for i in range(AtomNumber):
                fid.write("%12.6f %12.6f %12.6f\n"%(AtomPosition[0,i],AtomPosition[1,i],AtomPosition[2,i]))
            fid.close()
        os.system("cp POSCAR POSCAR.vasp")
        return

    def Intro_Screw_dipole(self):
        AtomNumber, H, Comm, AtomPosition = self.read_pos()
        xc1, xc2 = 0.25 * H[0,0], 0.75 * H[0,0] ; yc =  0.5 * H[1,1]
        for i in range(AtomNumber):
            dx1, dx2,  dy = AtomPosition[0, i] - xc1, AtomPosition[0, i] - xc2, AtomPosition[1, i] - yc
            theta1 = np.arctan2(dy, dx1)
            dz1 = self.screw_coeff * theta1
            theta2 = np.arctan2(dy, dx2)
            dz2 = self.screw_coeff * theta2
            AtomPosition[2, i] = AtomPosition[2, i] + dz1 - dz2
        with open("OUT", 'w') as fid:
            print >> fid , len(AtomPosition[0,:])
            print >> fid, len(AtomPosition[1,:])
            print >> fid , len(AtomPosition[2,:])
            fid.close()

        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.scatter(AtomPosition[0,:], AtomPosition[1,:], AtomPosition[2,:])
        plt.show()
        return

if __name__=="__main__":
    M = MD_ChangeBox()
    #M.Intro_Crack_xyz(1.4828,(-0.00809419040147+9.50392012258e-19j),(-0.00142961912234+3.90049109137e-19j),(2.00783401372e-20-0.00249498017278j),(9.49992985589e-19-0.00463795644713j),(1.11022302463e-16+1.74520621177j),(1.38777878078e-16+0.572998189701j))
    # M.Intro_Screw_dipole()
    #N = VA_ChangeBox()
    #N.volume_conserving_Ortho_strain(0.5)
    #N.volume_conserving_Mono_strain(0.5)
    #M.Intro_Screw_dipole()
    M.Mpi_Intro_Crack_xyz(1.4828,(-0.00809419040147+9.50392012258e-19j),(-0.00142961912234+3.90049109137e-19j),(2.00783401372e-20-0.00249498017278j),(9.49992985589e-19-0.00463795644713j),(1.11022302463e-16+1.74520621177j),(1.38777878078e-16+0.572998189701j))
