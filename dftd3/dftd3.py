#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

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
try:
    from .pars import *
except:
    from pars import *

#Python libararies
import random, sys, os, subprocess, string, math
from argparse import ArgumentParser
from cclib.io import ccread

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

FUNC_LIST = ["B1B95","B2GPPLYP","B3LYP","BHLYP","BLYP","BP86","BPBE","mPWLYP","PBE","PBE0","PW6B95","PWB6K","revPBE","TPSS","TPSS0","TPSSh","BOP","MPW1B95","MPWB1K","OLYP","OPBE","oTPSS","PBE38","PBEsol","REVSSB","SSB","B3PW91","BMK","CAMB3LYP","LCwPBE","M052X","M05","M062X","M06HF","M06L","M06","HCTH120","B2PLYP","DSDBLYP","TPSS","PWPB95","revPBE0","revPBE38","rPW86PBE"]
SUPPORTED_EXTENSIONS = set(('out', 'log', 'sdf', 'xyz', 'pdb', 'sdf'))

# Some useful arrays
periodictable = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
                 "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                 "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                 "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
                 "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
                 "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
                 "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
                 "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"]


def parseIntSet(nputstr=""):
    """
    Used to interpret a range of integers supplied at input (used to define fragment atoms)
    For example, the string '1-4, 7-9' will be expanded into a list of ints [1,2,3,4,7,8,9]
    """
    selection, invalid = set(), set()

    # tokens are comma seperated values
    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        if len(i) > 0:
            if i[:1] == "<":
                i = "1-%s"%(i[1:])
        try:
            # typically tokens are plain old integers
            selection.add(int(i))
        except:
            # if not, then it might be a range
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    first = token[0]
                    last = token[len(token)-1]
                    for x in range(first, last+1):
                        selection.add(x)
            except:
                # not an int and not a range...
                invalid.add(i)
    # Report invalid tokens before returning valid selection
    if len(invalid) > 0:
        print("Invalid set: " + str(invalid))
    return selection

## DFT derived values for diatomic cutoff radii from Grimme ##
## These are read from pars.py and converted from atomic units into Angstrom
r = [[0]*max_elem for x in range(max_elem)]

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
   c6ab = [[0]*max_elem for x in range(max_elem)]
   nlines = 32385

   for iat in range(0,max_elem):
      for jat in range(0,max_elem): c6ab[iat][jat]=[[0]*maxc for x in range(maxc)]
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

#verbose = None

## The computation of the D3 dispersion correction
class calcD3:
   def __init__(self, fileData, functional, damp='zero', s6=0.0, rs6=0.0, s8=0.0, a1=0.0, a2=0.0, abc=False, intermolecular=False, pairwise=False, verbose=False):

      atom_nums = fileData.atomnos.tolist()
      atomtype = [periodictable[atno] for atno in atom_nums]
      cartesians = fileData.atomcoords[-1].tolist()
      natom = len(atomtype)

      xco=[]; yco=[]; zco=[]
      for at in cartesians:
         xco.append(at[0]); yco.append(at[1]); zco.append(at[2])

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
      #print "\n           R0(Ang)        CN"
      #print "   #########################"
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
         #print "  ",(j+1), atomtype[j], "   %.5f" % dum, "   %.5f" % cn[j] #, C6jj, C8jj, C10jj
      #print "   #########################"

      icomp = [0]*100000; cc6ab = [0]*100000; r2ab = [0]*100000; dmp = [0]*100000

      ## Compute and output the individual components of the D3 energy correction ##
      #print "\n   Atoms  Types  C6            C8            E6              E8"
      if damp == "zero":
         if verbose: print("\n   D3-dispersion correction with zero-damping:", end=' ')
         if s6 == 0.0 or rs6 == 0.0 or s8 == 0.0:
            if functional != None:
               for parm in zero_parms:
                  if functional == parm[0]:
                     [s6,rs6,s8] = parm[1:4]
                     if verbose: print("detected", parm[0], "functional - using default zero-damping parameters")
            else:
                if verbose:
                    print("   WARNING: No functional information could be read!\n")
         else:
             if verbose: print(" manual parameters have been defined")
         if verbose: print("   Zero-damping parameters:", "s6 =",s6, "rs6 =", rs6, "s8 =",s8)

      if damp == "bj":
         if verbose: print("\n   D3-dispersion correction with Becke_Johnson damping:", end=' ')
         if s6 == 0.0 or s8 == 0.0 or a1 == 0.0 or a2 == 0.0:
            if functional != None:
               for parm in bj_parms:
                  if functional == parm[0]:
                     [s6,a1,s8,a2] = parm[1:5]
                     if verbose: print("detected", parm[0], "functional - using default BJ-damping parameters")
            else:
                if verbose: print("   WARNING: No functional information could be read!\n")
         else:
             if verbose: print(" manual parameters have been defined")
         if verbose: print("   BJ-damping parameters:", "s6 =",s6, "s8 =", s8, "a1 =",a1, "a2 =",a2)

      if verbose:
          if abc == False: print("   3-body term will not be calculated\n")
          else: print("   Including the Axilrod-Teller-Muto 3-body dispersion term\n")
          if intermolecular == True: print("   Only computing intermolecular dispersion interactions! This is not the total D3-correction\n")

      if intermolecular != False:
         mols = [0] * natom
         assign = [parseIntSet(im) for im in intermolecular.split(':')]
         for j, asn in enumerate(assign):
             for at in asn:
                 mols[int(at)-1] = j

      # loop over atom pairs
      for j in range(0,natom):
         for k in range(j+1,natom):
            scalefactor=1.0

            if intermolecular != False:
               if mols[j] == mols[k]:
                  scalefactor = 0
                  if verbose: print("   --- Ignoring interaction between atoms",(j+1), "and", (k+1))

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

               if pairwise == True and scalefactor != 0: print("   --- Pairwise interaction between atoms", (j+1), "and", (k+1),": Edisp =", "%.6f" % (self.attractive_r6_term + self.attractive_r8_term), "kcal/mol")

               self.attractive_r6_vdw = self.attractive_r6_vdw + self.attractive_r6_term
               self.attractive_r8_vdw = self.attractive_r8_vdw + self.attractive_r8_term

               jk=int(lin(k,j))
               icomp[jk] = 1
               cc6ab[jk] = math.sqrt(C6jk)
               r2ab[jk] = dist**2
               dmp[jk] = (1.0/rr)**(1.0/3.0)

      e63 = 0.0
      for iat in range(0,natom):
         for jat in range(0,natom):
            ij=int(lin(jat,iat))
            if icomp[ij]==1:
               for kat in range(jat,natom):
                  ik=int(lin(kat,iat))
                  jk=int(lin(kat,jat))

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

def main():
    parser = ArgumentParser()
    parser.add_argument("-v", dest="verbose", action="store_true", default=False, help="Turn on verbose printing")
    parser.add_argument("--damp", dest="damp", default="zero", type=str.lower, choices=('zero', 'bj'),
                        help="Type of D3-damping function (zero, bj)")
    parser.add_argument("--func", dest="functional", default=None,
                        help="Use default D3 parameters for this density functional")
    parser.add_argument("--s6", dest="s6", default=0.0, type=float, help="s6 parameter (used in zero and bj damping)")
    parser.add_argument("--rs6", dest="rs6", default=0.0, type=float, help="rs6 parameter used in zero damping")
    parser.add_argument("--s8", dest="s8", default=0.0, type=float, help="s8 parameter used in zero damping")
    parser.add_argument("--a1", dest="a1", default=0.0, type=float, help="a1 parameter used in bj damping")
    parser.add_argument("--kcal", dest="kcal", action="store_true", default=False, help="Print energies in kcal/mol")
    parser.add_argument("--a2", dest="a2", default=0.0, type=float, help="a2 parameter used in bj damping")
    parser.add_argument("--3body", dest="threebody", action="store_true", default=False, help="Turn on repulsive 3-body term")
    parser.add_argument("--pw", dest="pairwise", action="store_true", default=False,
                        help="Print dispersion terms between all interatomic pairs")
    parser.add_argument("--im", dest="intermolecular", type=str, default=False,
                        help="Compute only intermolecular dispersion terms")

    # Parse Arguments
    (options, args) = parser.parse_known_args()

    dft_functional = None
    if options.functional is not None:
        if options.functional in FUNC_LIST: dft_functional = options.functional
        else: print("\nUnable to match requested functional {} to stored parameters!\n".format(options.functional)); sys.exit()

    files = []
    for argv in sys.argv:
        if len(argv.split(".")) > 1:
            if argv.split(".")[-1] in SUPPORTED_EXTENSIONS:
                files.append(argv)
    if len(files)==0:
        print("\nNo valid files found!\n"); sys.exit()

    for file in files:
        try:
            data = ccread(file)

            if dft_functional == None:
                try:
                    if data.metadata['functional'] in FUNC_LIST: dft_functional = data.metadata['functional']
                except: pass

            fileD3 = calcD3(data, dft_functional, options.damp, options.s6, options.rs6, options.s8, options.a1, options.a2, options.threebody, options.intermolecular, options.pairwise, options.verbose)

            if not options.kcal:
                c6_term = fileD3.attractive_r6_vdw/autokcal
                c8_term = fileD3.attractive_r8_vdw/autokcal
                threebody_term = fileD3.repulsive_abc/autokcal

            else:
                c6_term = fileD3.attractive_r6_vdw
                c8_term = fileD3.attractive_r8_vdw
                threebody_term = fileD3.repulsive_abc

            # Output includes 3-body term
            if options.threebody == True:
                total_vdw = c6_term + c8_term + threebody_term
                if options.verbose: print('   {:<30} {:>13} {:>13} {:>13}'.format("Species", "D3(R6)", "D3(R8)", "3-body", "Total"))
                print('   {:<30} {:13.6f} {:13.6f} {:13.6f} {:13.6f}'.format(file, c6_term, c8_term, threebody_term, total_vdw))

            # Without 3-body term (default)
            else:
                total_vdw = c6_term + c8_term
                if options.verbose: print('   {:<30} {:>13} {:>13} {:>13}'.format("Species", "D3(R6)", " D3(R8)", "Total"))
                print('   {:<30} {:13.6f} {:13.6f} {:13.6f}'.format(file, c6_term, c8_term, total_vdw))
        except: pass
    print()

if __name__ == "__main__":
    main()
