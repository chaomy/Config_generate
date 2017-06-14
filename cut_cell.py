#!/usr/local/bin/python
# encoding: utf-8

import numpy as np
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f","--file",dest="filename")
parser.add_option("-s","--size", dest = "size", nargs = 3, type = "int",  help=" x y z size", default = "1")
(options,args) = parser.parse_args()

class cut(object):
    def __init__(self,filename,size):
        self.filename = filename
        self.size = size
    def geneposcar(self):
        BaseVectorx,BaseVectory,BaseVectorz = [],[],[]
        xx, yy, zz  = self.size[0], self.size[1], self.size[2]
        with open(self.filename) as fid :
            raw = fid.readlines()
        fid.close()
        #------------read file --------------
        AtomName = raw[0]
        lattice_constant = raw[1]
        for i in range(2,5):
            BaseVectorx.append(float(raw[i].split()[i-2])/self.size[0])
            BaseVectory.append(float(raw[i].split()[i-2])/self.size[1])
            BaseVectorz.append(float(raw[i].split()[i-2])/self.size[2])
        Element = raw[5]
        ElementNum = []
        for i in raw[6].split():
            ElementNum.append(int(i))
        xx, yy, zz  = [],[],[]
        for line in raw[8:-1]:
            xx.append(float(line.split()[0])) ; yy.append(float(line.split()[1])) ; zz.append(float(line.split()[2]))
        #-----------cut -----------------
        Xc, Yc, Zc = 1./self.size[0], 1./self.size[1], 1./self.size[2]
        xnew , ynew, znew = [],[],[]
        for i in range(len(xx)):
            if xx[i] < Xc and yy[i] < Yc and zz[i] < Zc :
                xnew.append(xx[i] * self.size[0]);  ynew.append(yy[i] * self.size[1]) ;  znew.append(zz[i]* self.size[2])
        ElementNumnew = []
        Coeff = int (self.size[0] * self.size[1] * self.size[2])
        for i in ElementNum:
            ElementNumnew.append(i/Coeff)
        #----------output ---------------
        with open("POSCAR", 'w') as fid:
            fid.write("%s"%(AtomName))
            fid.write("%s"%(lattice_constant))
            fid.write("%7.5f\t%7.5f\t%7.5f\n"%(BaseVectorx[0],BaseVectorx[1],BaseVectorx[2]))
            fid.write("%7.5f\t%7.5f\t%7.5f\n"%(BaseVectory[0],BaseVectory[1],BaseVectory[2]))
            fid.write("%7.5f\t%7.5f\t%7.5f\n"%(BaseVectorz[0],BaseVectorz[1],BaseVectorz[2]))
            fid.write("%s"%(Element))
            fid.write("%d\t%d\t%d\n"%(ElementNumnew[0],ElementNumnew[1],ElementNumnew[2]))
            fid.write("Direct\n")
            for i in range(len(xnew)):
                fid.write("%7.5f\t%7.5f\t%7.5f\n"%(xnew[i],ynew[i],znew[i]))
M=cut(options.filename,options.size)
M.geneposcar()
