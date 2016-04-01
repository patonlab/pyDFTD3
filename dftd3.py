#!/usr/bin/python

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#######################################################################
#                          dftd3.py                                   #
#                                                                     #
#    This was a little exercise to translate Grimme's D3              #
#    Fortran code into Python. There is no new science as such!       #
#    It is possible to implement D3 corrections directly within most  #
#    electronic structure packages so the only possible use for this  #
#    is pedagogic or to look at individual terms of the total         #
#    D3-correction between a pair of atoms or molecules.              #
#    This code will read a Gaussian formatted input/output file and   #
#    compute the D3-density independent dispersion terms without      #
#    modification. Zero and Becke-Johnson damping schemes are both    #
#    implemented.                                                     #
#######################################################################
#######  Written by:  Rob Paton and Kelvin Jackson ####################
#######  Last modified:  Mar 20, 2016 #################################
#######################################################################

# Dependent on parameter file
from pars import *

# For reading Gaussian formatted input/output files
from ccParse import *

#Python libararies
import random, sys, os, commands, string, math

## Check for integer when parsing ##
def is_number(s):
    try: int(s); return True
    except ValueError: return False

## Arrays for attractive and repulsive interactions ##
attractive_vdw=[0]
repulsive_vdw=[0]
total_vdw=[0]

## Functional Specific D3 parameters
rs8 = 1.0
repfac= 1.0

## Distance and energy conversion factors ##
autoang = 0.52917726
autokcal = 627.509541
c6conv=(0.001/2625.4999)/(0.052917726**6)

## Global D3 parameters ##
## Exponents used in distance dependent damping factors for R6, R8 and R10 terms
alpha6 = 14
alpha8 = alpha6 + 2
alpha10 = alpha8 + 2

## Constants used to determine fractional connectivities between 2 atoms:
## k1 is the exponent used in summation, k2 is used a fraction of the summed single-bond radii
k1 = 16.0
k2 = 4.0/3.0
k3 = -4.0

## D3 is parameterized up to element 94
max_elem = 94
## maximum connectivity
maxc = 5

## From connectivity, establish if there is more than one molecule
def getMollist(bondmatrix,startatom):
   # The list of atoms in a molecule
   atomlist=[]
   atomlist.append(startatom)
   molecule1=[]
   nextlot=[]
   count = 0

   while count<100:
      nextlot=[]
      for atom in atomlist:
         for i in range(0,len(bondmatrix[atom])):
            if bondmatrix[atom][i] == 1:
               alreadyfound = 0
               for at in atomlist:
                  if i == at: alreadyfound = 1
               if alreadyfound == 0: atomlist.append(i)
      count=count+1

   return atomlist

## DFT derived values for diatomic cutoff radii from Grimme ##
## These are read from pars.py and converted from atomic units into Angstrom
r = [[0]*max_elem for x in xrange(max_elem)]

k=0
for i in range(0,max_elem):
   for j in range(0,i+1):
      r[i][j]=r0ab[k]/autoang
      r[j][i]=r0ab[k]/autoang
      k=k+1

## PBE0/def2-QZVP atomic values for multipole coefficients read from pars.py ##
for i in range(0,max_elem):
   dum=0.5*r2r4[i]*float(i+1)**0.5
   r2r4[i]=math.pow(dum,0.5)

## Reference systems are read in to compute coordination number dependent dispersion coefficients
def copyc6(max_elem, maxc):
   c6ab = [[0]*max_elem for x in xrange(max_elem)]
   nlines = 32385

   for iat in range(0,max_elem):
      for jat in range(0,max_elem): c6ab[iat][jat]=[[0]*maxc for x in xrange(maxc)]
   kk=0

   for nn in range(0,nlines):
      kk=(nn*5)
      iadr=0
      jadr=0
      iat=int(pars[kk+1])-1
      jat=int(pars[kk+2])-1

      while iat > 99:
         iadr=iadr+1
         iat=iat-100
      while jat > 99:
         jadr=jadr+1
         jat=jat-100

      c6ab[iat][jat][iadr][jadr]=[]
      c6ab[iat][jat][iadr][jadr].append(pars[kk])
      c6ab[iat][jat][iadr][jadr].append(pars[kk+3])
      c6ab[iat][jat][iadr][jadr].append(pars[kk+4])

      c6ab[jat][iat][jadr][iadr]=[]
      c6ab[jat][iat][jadr][iadr].append(pars[kk])
      c6ab[jat][iat][jadr][iadr].append(pars[kk+4])
      c6ab[jat][iat][jadr][iadr].append(pars[kk+3])

   return c6ab

# Obtain the C6 coefficient for the interaction between atoms a and b
def getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,a,b):
   for i in range(0,max_elem):
      if atomtype[a].find(elements[i])>-1:iat=i
      if atomtype[b].find(elements[i])>-1:jat=i

   c6mem = -1.0E99
   rsum = 0.0
   csum = 0.0
   c6  = 0.0

   for i in range(0,mxc[iat]):
      for j in range(0,mxc[jat]):
         if isinstance(c6ab[iat][jat][i][j], (list, tuple)):
            c6=c6ab[iat][jat][i][j][0]
            if c6>0:
               c6mem=c6
               cn1=c6ab[iat][jat][i][j][1]
               cn2=c6ab[iat][jat][i][j][2]

               r=(cn1-cn[a])**2+(cn2-cn[b])**2
               tmp1=math.exp(k3*r)
               rsum=rsum+tmp1
               csum=csum+tmp1*c6

   if(rsum>0): c6=csum/rsum
   else: c6=c6mem

   return c6

# Calculation of atomic coordination numbers
def ncoord(natom, rcov, atomtype, xco, yco, zco, max_elem, autoang, k1, k2):
   cn =[]
   for i in range(0,natom):
      xn = 0.0
      for iat in range(0,natom):
         if iat != i:
            dx = xco[iat] - xco[i]
            dy = yco[iat] - yco[i]
            dz = zco[iat] - zco[i]
            r2 = dx*dx+dy*dy+dz*dz
            r = math.pow(r2,0.5)
            r = r

            for k in range(0,max_elem):
               if atomtype[i].find(elements[k])>-1:Zi=k
               if atomtype[iat].find(elements[k])>-1:Ziat=k

            rco = rcov[Zi]+rcov[Ziat]
            rco = rco*k2
            rr=rco/r
            damp=1.0/(1.0+math.exp(-k1*(rr-1.0)))
            xn=xn+damp
      cn.append(xn)

   return cn

# linear interpolation
def lin(i1,i2):
   idum1=max(i1,i2)
   idum2=min(i1,i2)
   lin=idum2+idum1*(idum1-1)/2
   return lin

## Get from pars.py
c6ab = copyc6(max_elem, maxc)

## The periodic table...
elm=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

## Published parameters (S6, RS6, S8) for zero-damping optimized for different functionals
zero_parms = [["B2PLYP",	0.6400,	1.4270,	1.0220],["B97",	1.0000,	0.8920,	0.9090],["B3LYP",	1.0000,	1.2610,	1.7030],["BLYP",	1.0000,	1.0940,	1.6820],["PBE1PBE",	1.0000,	1.2870,	0.9280],["TPSSTPSS",	1.0000,	1.1660,	1.1050],["PBEPBE",	1.0000,	1.2170,	0.7220],["BP86",	1.0000,	1.1390,	1.6830],["BPBE",	1.0000,	1.0870,	2.0330],["B3PW91",	1.0000,	1.1760,	1.7750],["BMK",	1.0000,	1.9310,	2.1680],["CAM-B3LYP",	1.0000,	1.3780,	1.2170],["LC-wPBE",	1.0000,	1.3550,	1.2790],["M05",	1.0000,	1.3730,	0.5950],["M052X",	1.0000,	1.4170,	0.0000],["M06L",	1.0000,	1.5810,	0.0000],["M06",	1.0000,	1.3250,	0.0000],["M062X",	1.0000,	1.6190,	0.0000],["M06HF",	1.0000,	1.4460,	0.0000]]

## Published parameters (S6, S8, a1, a2) for Becke-Johnson-damping optimized for different functionals
bj_parms = [["B2PLYP",	0.6400,	0.9147,	0.3065,	5.0570],["B97",	1.0000,	2.2609,	0.5545,	3.2297],["B3LYP",	1.0000,	1.9889,	0.3981,	4.4211],["BLYP",	1.0000,	2.6996,	0.4298,	4.2359],["PBE1PBE",	1.0000,	1.2177,	0.4145,	4.8593],["TPSSTPSS",	1.0000,	1.9435,	0.4535,	4.4752],["PBEPBE",	1.0000,	0.7875,	0.4289,	4.4407],["BP86",	1.0000,	3.2822,	0.3946,	4.8516],["BPBE",	1.0000,	4.0728,	0.4567,	4.3908],["B3PW91",	1.0000,	2.8524,	0.4312,	4.4693],["BMK",	1.0000,	2.0860,	0.1940,	5.9197],["CAM-B3LYP",	1.0000,	2.0674,	0.3708,	5.4743],["LC-wPBE",	1.0000,	1.8541,	0.3919,	5.0897]]

## The computation of the D3 dispersion correction
class calcD3:
   def __init__(self, file, s6, rs6, s8, a1, a2, damp, abc, intermolecular, pairwise):

      ## Use ccParse to get the Cartesian coordinates from Gaussian input/output files
      if len(file.split(".com"))>1 or len(file.split(".gjf"))>1: fileData = getinData(file)
      if len(file.split(".out"))>1 or len(file.split(".log"))>1: fileData = getoutData(file)

      ## Arrays for atoms and Cartesian coordinates ##
      atomtype = fileData.ATOMTYPES
      natom = len(atomtype)

      xco=[]; yco=[]; zco=[]
      for at in fileData.CARTESIANS:
         xco.append(at[0]); yco.append(at[1]); zco.append(at[2])

      ## In case something clever needs to be done wrt inter and intramolecular interactions
      if hasattr(fileData,"BONDINDEX"):
         molAatoms = getMollist(fileData.BONDINDEX,0)
         mols = []
         for j in range(0,natom):
            mols.append(0)
            for atom in molAatoms:
               if atom == j: mols[j] = 1

      ## Names are pretty obvious...
      self.attractive_r6_vdw = 0.0; self.attractive_r8_vdw = 0.0; self.repulsive_abc = 0.0

      mxc=[0]
      for j in range(0,max_elem):
         mxc.append(0)
         for k in range(0,natom):
            if atomtype[k].find(elements[j])>-1:
               for l in range(0,maxc):
                  if isinstance(c6ab[j][j][l][l], (list, tuple)):
                     if c6ab[j][j][l][l][0]>0: mxc[j]=mxc[j]+1
               break

      ## Coordination number based on covalent radii
      cn = ncoord(natom, rcov, atomtype, xco, yco, zco, max_elem, autoang, k1, k2)

      ## C6 - Need to calculate these from fractional coordination
      print "\n           R0(Ang)        CN"
      print "   #########################"
      x=0
      for j in range(0,natom):
         C6jj = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,j)

         for k in range(0,natom):
            dum = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,k)
            x=x+dum

         for k in range(0,max_elem):
            if atomtype[j].find(elements[k])>-1:z=k

         dum = 0.5*autoang*r[z][z]

         C8jj = 3.0*C6jj*math.pow(r2r4[z],2.0)
         C10jj=49.0/40.0 * math.pow(C8jj,2.0)/C6jj
         print "  ",(j+1), atomtype[j], "   %.5f" % dum, "   %.5f" % cn[j] #, C6jj, C8jj, C10jj
      print "   #########################"

      icomp = [0]*100000; cc6ab = [0]*100000; r2ab = [0]*100000; dmp = [0]*100000
                            
      ## Compute and output the individual components of the D3 energy correction ##
      #print "\n   Atoms  Types  C6            C8            E6              E8"
      if damp == "zero":
         print "\n   D3-dispersion correction with zero-damping:",
         if s6 == 0.0 or rs6 == 0.0 or s8 == 0.0:
            if fileData.FUNCTIONAL != "none":
               for parm in zero_parms:
                  if fileData.FUNCTIONAL == parm[0]:
                     [s6,rs6,s8] = parm[1:4]
                     print "detected", parm[0], "functional - using default zero-damping parameters"
            else: print "   WARNING: Damping parameters not specified and no functional could be read!\n"; sys.exit()
         else: print " manual parameters have been defined"
         print "   Zero-damping parameters:", "s6 =",s6, "rs6 =", rs6, "s8 =",s8
     
      if damp == "bj":
         print "\n   D3-dispersion correction with Becke_Johnson damping:",
         if s6 == 0.0 or s8 == 0.0 or a1 == 0.0 or a2 == 0.0:
            if fileData.FUNCTIONAL != "none":
               for parm in bj_parms:
                  if fileData.FUNCTIONAL == parm[0]:
                     [s6,s8,a1,a2] = parm[1:5]
                     print "detected", parm[0], "functional - using default BJ-damping parameters"
            else: print "   WARNING: Damping parameters not specified and no functional could be read!\n"; sys.exit()
         else: print " manual parameters have been defined"
         print "   BJ-damping parameters:", "s6 =",s6, "s8 =", s8, "a1 =",a1, "a2 =",a2

      if abc == "off": print "   3-body term will not be calculated\n"
      else: print "   Including the Axilrod-Teller-Muto 3-body dispersion term\n"
      if intermolecular == "on": print "   Only computing intermolecular dispersion interactions! This is not the total D3-correction\n"

      for j in range(0,natom):
         ## This could be used to 'switch off' dispersion between bonded or geminal atoms ##
         scaling = "off"
         for k in range(j+1,natom):
            scalefactor=1.0

            if intermolecular == "on":
               if mols[j] == mols[k]:
                  scalefactor = 0
                  print "   --- Ignoring interaction between atoms",(j+1), "and", (k+1)

            if scaling=="on" and hasattr(fileData,"BONDINDEX"):
               if fileData.BONDINDEX[j][k]==1: scalefactor = 0
               for l in range (0,natom):
                  if fileData.BONDINDEX[j][l] != 0 and fileData.BONDINDEX[k][l]!=0 and j!=k and fileData.BONDINDEX[j][k]==0: scalefactor = 0
                  for m in range (0,natom):
                          if fileData.BONDINDEX[j][l] != 0 and fileData.BONDINDEX[l][m]!=0 and fileData.BONDINDEX[k][m]!=0 and j!=m and k!=l and fileData.BONDINDEX[j][m]==0: scalefactor=1/1.2

            if k>j:
               ## Pythagoras in 3D to work out distance ##
               xdist = xco[j]-xco[k]; ydist = yco[j]-yco[k]; zdist = zco[j]-zco[k]
               totdist = math.pow(xdist,2)+math.pow(ydist,2)+math.pow(zdist,2)
               totdist=math.sqrt(totdist)

               C6jk = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,k)

               ## C8 parameters depend on C6 recursively
               for l in range(0,max_elem):
                  if atomtype[j].find(elements[l])>-1:atomA=l
                  if atomtype[k].find(elements[l])>-1:atomB=l

               C8jk = 3.0*C6jk*r2r4[atomA]*r2r4[atomB]
               C10jk=49.0/40.0 * math.pow(C8jk,2.0)/C6jk

               # Evaluation of the attractive term dependent on R^-6 and R^-8
               if damp == "zero":
                  dist=totdist/autoang
                  rr = r[atomA][atomB]/dist
                  tmp1 = rs6*rr
                  damp6 = 1/(1+6*math.pow(tmp1,alpha6))
                  tmp2 = rs8*rr
                  damp8 = 1/(1+6*math.pow(tmp2,alpha8))

                  self.attractive_r6_term = -s6*C6jk*damp6/math.pow(dist,6)*autokcal*scalefactor
                  self.attractive_r8_term = -s8*C8jk*damp8/math.pow(dist,8)*autokcal*scalefactor
               
               if damp == "bj":
                  dist=totdist/autoang
                  rr = r[atomA][atomB]/dist
                  rr = math.pow((C8jk / C6jk), 0.5)
                  tmp1 = a1 * rr + a2
                  damp6 = math.pow(tmp1,6)
                  damp8 = math.pow(tmp1,8)
                  
                  self.attractive_r6_term = -s6*C6jk/(math.pow(dist,6)+damp6)*autokcal*scalefactor
                  self.attractive_r8_term = -s8*C8jk/(math.pow(dist,8)+damp8)*autokcal*scalefactor
                  
               if pairwise == "on" and scalefactor != 0: print "   --- Pairwise interaction between atoms", (j+1), "and", (k+1),": Edisp =", "%.6f" % (self.attractive_r6_term + self.attractive_r8_term), "kcal/mol"

               self.attractive_r6_vdw = self.attractive_r6_vdw + self.attractive_r6_term
               self.attractive_r8_vdw = self.attractive_r8_vdw + self.attractive_r8_term

               jk=lin(k,j)
               icomp[jk] = 1
               cc6ab[jk] = math.sqrt(C6jk)
               r2ab[jk] = dist**2
               dmp[jk] = (1.0/rr)**(1.0/3.0)

      e63 = 0.0
      for iat in range(0,natom):
         for jat in range(0,natom):
            ij=lin(jat,iat)
            if icomp[ij]==1:
               for kat in range(jat,natom):
                  ik=lin(kat,iat)
                  jk=lin(kat,jat)

                  if kat>jat and jat>iat and icomp[ik] != 0 and icomp[jk] != 0:
                     rav=(4.0/3.0)/(dmp[ik]*dmp[jk]*dmp[ij])
                     tmp=1.0/( 1.0+6.0*rav**alpha6 )

                     c9=cc6ab[ij]*cc6ab[ik]*cc6ab[jk]
                     d2 = [0]*3
                     d2[0]=r2ab[ij]
                     d2[1]=r2ab[jk]
                     d2[2]=r2ab[ik]
                     t1 = (d2[0]+d2[1]-d2[2])/math.sqrt(d2[0]*d2[1])
                     t2 = (d2[0]+d2[2]-d2[1])/math.sqrt(d2[0]*d2[2])
                     t3 = (d2[2]+d2[1]-d2[0])/math.sqrt(d2[1]*d2[2])
                     ang=0.375*t1*t2*t3+1.0
                     e63=e63+tmp*c9*ang/(d2[0]*d2[1]*d2[2])**1.50

      self.repulsive_abc_term = s6 * e63 * autokcal
      self.repulsive_abc = self.repulsive_abc + self.repulsive_abc_term

if __name__ == "__main__":
   # Takes arguments: (1) damping style, (2) s6, (3) rs6, (4) s8, (5) 3-body on/off, (6) input file(s)
   files = []
   damp = "zero"; s6 = 0.0; rs6 = 0.0; s8 = 0.0; bj_a1 = 0.0; bj_a2 = 0.0; abc_term = "off"; intermolecular = "off"; pairwise = "off"
   if len(sys.argv) > 1:
      for i in range(1,len(sys.argv)):
         if sys.argv[i] == "-damp": damp = (sys.argv[i+1])
         elif sys.argv[i] == "-s6": s6 = float(sys.argv[i+1])
         elif sys.argv[i] == "-rs6": rs6 = float(sys.argv[i+1])
         elif sys.argv[i] == "-s8": s8 = float(sys.argv[i+1])
         elif sys.argv[i] == "-a1": bj_a1 = float(sys.argv[i+1])
         elif sys.argv[i] == "-a2": bj_a2 = float(sys.argv[i+1])
         elif sys.argv[i] == "-3body": abc_term = (sys.argv[i+1])
         elif sys.argv[i] == "-im": intermolecular = (sys.argv[i+1])
         elif sys.argv[i] == "-pw": pairwise = (sys.argv[i+1])
         else:
            if len(sys.argv[i].split(".")) > 1:
               if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log" or sys.argv[i].split(".")[1] == "com" or sys.argv[i].split(".")[1] == "gjf":
                  files.append(sys.argv[i])

   else: print "\nWrong number of arguments used. Correct format: dftd3.py (-damp zero/bj) (-s6 val) (-rs6 val) (-s8 val) (-a1 val) (-a2 val) (-im on/off) (-pw on/off)file(s)\n"; sys.exit()

   for file in files:
      fileD3 = calcD3(file, s6, rs6, s8, bj_a1, bj_a2, damp, abc_term, intermolecular, pairwise)
      attractive_r6_vdw = fileD3.attractive_r6_vdw/autokcal
      attractive_r8_vdw = fileD3.attractive_r8_vdw/autokcal

      # Output includes 3-body term
      if abc_term == "on":
         repulsive_abc = fileD3.repulsive_abc/autokcal
         total_vdw = attractive_r6_vdw + attractive_r8_vdw + repulsive_abc
         print "\n", " ".rjust(30), "    D3(R6)".rjust(12), "    D3(R8)".rjust(12),"    D3(3-body)".rjust(12), "    Total (au)".rjust(12)
         print "  ",file.ljust(30), "   %.8f" % attractive_r6_vdw, "   %.8f" % attractive_r8_vdw,"   %.8f" % repulsive_abc, "   %.8f" % total_vdw, "\n"

      # Without 3-body term (default)
      else:
         total_vdw = attractive_r6_vdw + attractive_r8_vdw
         print "\n", " ".rjust(30), "    D3(R6)".rjust(12), "    D3(R8)".rjust(12), "    Total (au)".rjust(12)
         print "  ",file.ljust(30), "   %.8f" % attractive_r6_vdw, "   %.8f" % attractive_r8_vdw, "   %.8f" % total_vdw, "\n"

