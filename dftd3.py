#!/usr/bin/python

###                   ###                     ###          ###      ###
###                   ###                     ###          ###      ###
    #####b.   ####b.  ###### .d##b.  #####b.  ###  ####b.  #####b.  ###
### ### "##b     "##b ###   d##""##b ### "##b ###     "##b ### "##b ###
### ###  ### .d###### ###   ###  ### ###  ### ### .d###### ###  ###
### ### d##P ###  ### Y##b. Y##..##P ###  ### ### ###  ### ### d##P ###
### #####P"  "Y######  "Y### "Y##P"  ###  ### ### "Y###### #####P"  ###
    ###                                                             
    ###

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
#                        eval_D3.py                                   #
#     Kelvin Jackson & Robert Paton, University of Oxford, 2012       #
#                                                                     #
#    This is an attempt to translate the original D3 Fortran code of  #
#    Grimme, as described in Stefan Grimme, Jens Antony, Stephan      #
#    Ehrlich, and Helge Krieg, J. Chem. Phys. 132, 154104 (2010)      #
#                                                                     #
#    For a Gaussian formatted input/output file this program will     #
#    compute Grimme's D3 attractive interaction                       #
#    All attractive potentials, cutoffs, damping, connectivity are    #
#    taken directly from Grimme's DFT-D3 program                      #
#    Where connectivity is available (*com), it can be used to scale  #
#   interactions and compute steric repulsions (under devlopment)     #
#######################################################################
#######  Written by:  Rob Paton #######################################
#######  Last modified:  Mar 20, 2013 #################################
#######################################################################


#  To be added:  Becke-Johnson Damping

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

## Conversion factors ##
autoang = 0.52917726
autokcal = 627.509541
c6conv=(0.001/2625.4999)/(0.052917726**6)

## Global D3 parameters ##
## Exponents used in distance dependent damping factors for R6, R8 and R10 terms
alpha6 = 14
alpha8 = alpha6 + 2
alpha10 = alpha8 + 2
## Constants used to determine the fractional connectivity between 2 atoms:
## k1 is the exponent used in summation, k2 is used a fraction of the summed single-bond radii
k1 = 16.0
k2 = 4.0/3.0
k3 = -4.0
## D3 is parameterized up to element 94
max_elem = 94
## maximum connectivity
maxc = 5

## Work out the atoms in the same molecules
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
			#print atom, bondmatrix[atom]
			
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
		for jat in range(0,max_elem):
			c6ab[iat][jat]=[[0]*maxc for x in xrange(maxc)]
	
	kk=0
	for nn in range(0,nlines):
		kk=(nn*5)
		iadr=0
		jadr=0
		#print pars[kk], pars[kk+1], pars[kk+2]
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
		
	if(rsum>0):
		c6=csum/rsum
	else:
		c6=c6mem	
	return c6

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
				#print Zi, Ziat,r, rcov[Zi], rcov[Ziat], rco,rr, damp
				xn=xn+damp
		
		cn.append(xn)
	return cn

def lin(i1,i2):
	idum1=max(i1,i2)
	idum2=min(i1,i2)
	lin=idum2+idum1*(idum1-1)/2
	return lin

## Get from pars.py
c6ab = copyc6(max_elem, maxc)

## The periodic table...
elm=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

## The computation of the D3 dispersion correction
class calcD3:
	def __init__(self, file, s6, rs6, s8):
		
		## Use ccParse to get the geomtetric data from file
		if len(file.split(".com"))>1:
			file = file.split(".com")[0]
			fileData = getinData(file)
				
		if len(file.split(".out"))>1 or len(file.split(".log"))>1:
			file = file.split(".")[0]
			fileData = getoutData(file)
			
			if hasattr(fileData,"FUNCTIONAL"):
				#if fileData.FUNCTIONAL == "B3LYP": s6 = 1.0000; rs6 = 1.2610; s8 = 1.7030; #print "   \no  Using default B3LYP D3 parameters:",
				if fileData.FUNCTIONAL == "BP86": s6 = 1.0000; rs6 = 1.139; s8 = 1.683; print "   \no  Using default BP86 D3 parameters:",
				if fileData.FUNCTIONAL == "B2PLYP": s6 = 0.5000; rs6 = 1.551; s8 = 1.1090; print "   \no  Using default B2PLYP/TZVPP D3 parameters:",
				if fileData.FUNCTIONAL == "M06-2X": s6 = 1.0000; rs6 = 1.6190; s8 = 0.0000; print "   \no  Using default M06-2X D3 parameters:",
				if fileData.FUNCTIONAL == "M06L": s6 = 1.0000; rs6 = 1.5810; s8 = 0.0000; print "   \no  Using default M06L D3 parameters:",
				if fileData.FUNCTIONAL == "M06": s6 = 1.0000; rs6 = 1.3250; s8 = 0.0000; print "   \no  Using default M06 D3 parameters:",
				if fileData.FUNCTIONAL == "B97D": s6 = 1.0000; rs6 = 0.8920; s8 = 0.9090; print "   \no  Using default B97-D D3 parameters:",
	
		#print "s6 =",s6, "rs6 = ", rs6, "s8 =",s8
		## Arrays for atoms and Cartesian coordinates ##
		atomtype = fileData.ATOMTYPES
		natom = len(atomtype)
		
		xco=[]; yco=[]; zco=[]
		for at in fileData.CARTESIANS:
			xco.append(at[0])
			yco.append(at[1])
			zco.append(at[2])

		## In case something clever needs to be done wrt inter and intramolecular interaction
		if hasattr(fileData,"BONDINDEX"):
			molAatoms = getMollist(fileData.BONDINDEX,0)
			mols = []
			for j in range(0,natom):
				mols.append(0)
				for atom in molAatoms:
					if atom == j: mols[j] = 1
	
		## Names are pretty obvious...
		self.repulsive_vdw = 0.0; self.attractive_r6_vdw = 0.0; self.attractive_r8_vdw = 0.0; self.repulsive_abc = 0.0
		
		#print "o  Using the following C6 cooefficients"
		#print "   ID   Z    CN        C6"
		mxc=[0]
		for j in range(0,max_elem):
			mxc.append(0)
			for k in range(0,natom):
				if atomtype[k].find(elements[j])>-1:
					for l in range(0,maxc):
						if isinstance(c6ab[j][j][l][l], (list, tuple)):
							if c6ab[j][j][l][l][0]>0:
								#print "  ", atomtype[k],"  ", (j+1),"  ", c6ab[j][j][l][l][1],"  ", c6ab[j][j][l][l][0]
								mxc[j]=mxc[j]+1
					break

		## Coordination number based on covalent radii
		cn = ncoord(natom, rcov, atomtype, xco, yco, zco, max_elem, autoang, k1, k2)

		## C6 - Need to calculate these from fractional coordination
		print "\n   #                XYZ [au]                   R0(AA) [Ang.]  CN          C6(AA)     C8(AA)   C10(AA) [au]"

		x=0
		for j in range(0,natom):
			
			C6jj = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,j)
			
			for k in range(0,natom):
				dum = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,k)
				x=x+dum
			#print j, k, x
			
			for k in range(0,max_elem):
				if atomtype[j].find(elements[k])>-1:z=k
			
			dum = 0.5*autoang*r[z][z]
			
			C8jj = 3.0*C6jj*math.pow(r2r4[z],2.0)
			C10jj=49.0/40.0 * math.pow(C8jj,2.0)/C6jj
			print "  ",(j+1), xco[j], yco[j], zco[j], atomtype[j], dum, cn[j], C6jj, C8jj, C10jj
		
		#print "\n   Molecular C6(AA) [au] =   ", x

		icomp = [0]*100000; cc6ab = [0]*100000; r2ab = [0]*100000; dmp = [0]*100000
					
		## Compute and output the individual components of the D3 energy correction ##
		#print "\n   Atoms  Types  C6            C8            E6              E8"
		for j in range(0,natom):
			## This could be used to 'switch off' dispersion between bonded or geminal atoms ##
			for k in range(j+1,natom):
				scalefactor=1.0
				vdw=0
				
				if hasattr(fileData,"BONDINDEX"):
					if fileData.BONDINDEX[j][k]==1: vdw=0; scalefactor = 0
					for l in range (0,natom):
						if fileData.BONDINDEX[j][l] != 0 and fileData.BONDINDEX[k][l]!=0 and j!=k and fileData.BONDINDEX[j][k]==0: vdw=0; scalefactor = 0
						for m in range (0,natom):
							if fileData.BONDINDEX[j][l] != 0 and fileData.BONDINDEX[l][m]!=0 and fileData.BONDINDEX[k][m]!=0 and j!=m and k!=l and fileData.BONDINDEX[j][m]==0: scalefactor=1/1.2
				
				
				if k>j:
					#print "  ", (j+1), (k+1),"  ", atomtype[j], atomtype[k],"  ", scalefactor
					
					## Pythagoras in 3D to work out distance ##
					xdist = xco[j]-xco[k]
					ydist = yco[j]-yco[k]
					zdist = zco[j]-zco[k]
					totdist = math.pow(xdist,2)+math.pow(ydist,2)+math.pow(zdist,2)
					totdist=math.sqrt(totdist)
					
					
					C6jk = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,k)
					
					## C8 parameters depend on C6 recursively
					for l in range(0,max_elem):
						if atomtype[j].find(elements[l])>-1:atomA=l
						if atomtype[k].find(elements[l])>-1:atomB=l
					
					C8jk = 3.0*C6jk*r2r4[atomA]*r2r4[atomB]
					C10jk=49.0/40.0 * math.pow(C8jk,2.0)/C6jk
					
					#print C6jk, C8jk,
					
					dist=totdist/autoang
					rr = r[atomA][atomB]/dist
					tmp1 = rs6*rr
					damp6 = 1/(1+6*math.pow(tmp1,alpha6))
					#print dist, r[atomA][atomB],rr, tmp1, damp6,
					tmp2 = rs8*rr
					damp8 = 1/(1+6*math.pow(tmp2,alpha8))
						
					## Repulsive coefficient
					R6jk = 4.5/math.pow(rs6,12)
					
					# Example for a repulsive potential dependent on R^-12
					# If two atoms are bonded then this term is zero
					if vdw == 1 and len(fileData.BONDINDEX) > 1:
						self.repulsive_vdw_term = R6jk*(1-damp6)/math.pow(1/tmp1,12)*math.pow(repfac,12)*scalefactor
					else: self.repulsive_vdw_term = 0.0
					
					#print R6jk, 1/math.pow(1/tmp1,12), math.pow(repfac,12)
					
					# Evaluation of the attractive term dependent on R^-6
					self.attractive_r6_term = -s6*C6jk*damp6/math.pow(dist,6)*autokcal*scalefactor
					
					# Evaluation of the attractive term dependent on R^-8
					self.attractive_r8_term = -s8*C8jk*damp8/math.pow(dist,8)*autokcal*scalefactor
					#self.attractive_r8_term = 0.0

					self.repulsive_vdw = self.repulsive_vdw + self.repulsive_vdw_term
					self.attractive_r6_vdw = self.attractive_r6_vdw + self.attractive_r6_term
					self.attractive_r8_vdw = self.attractive_r8_vdw + self.attractive_r8_term
					#print self.attractive_r6_term, self.attractive_r8_term, self.repulsive_vdw_term

	
					jk=lin(k,j)
					icomp[jk] = 1
					cc6ab[jk] = math.sqrt(C6jk)
					r2ab[jk] = dist**2
					dmp[jk] = (1.0/rr)**(1.0/3.0)
					#print j, k, jk

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
	
	# Takes arguments: (1) s6, (2) rs6, (3) s8, (4) input file(s)
	abc = 1
	files = []
	if len(sys.argv) > 4:
		s6 = float(sys.argv[1])
		rs6 = float(sys.argv[2])
		s8 = float(sys.argv[3])
		for i in range(4,len(sys.argv)):
			files.append(sys.argv[i])
	else:
		print "\nWrong number of arguments used. Correct format: eval_D3 s6 rs6 s8 file(s)\n"
		sys.exit()

	for file in files:
		abc = 0
		fileD3 = calcD3(file, s6, rs6, s8)
		attractive_r6_vdw = fileD3.attractive_r6_vdw/autokcal
		attractive_r8_vdw = fileD3.attractive_r8_vdw/autokcal
		if abc == 1: repulsive_abc = fileD3.repulsive_abc/autokcal
		else: repulsive_abc = 0.0
		total_vdw = attractive_r6_vdw + attractive_r8_vdw + repulsive_abc
		#print "   Breakdown   Attractive-R6   Attractive-R8   Repulsive-3-Body   Total   (Hartree)"
		print "  ",file, attractive_r6_vdw, attractive_r8_vdw,repulsive_abc, total_vdw
