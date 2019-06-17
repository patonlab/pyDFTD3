pyDFTD3
======

This program will compute the D3-dispersion energy, as developed by Grimme, for one or a series of Gaussian-formatted input or output file(s) in Cartesian coordinates. The dispersion energies will be identical to those produced by the original Fortran version of dftd3, or from the Gaussian program itself, provided the same damping parameters are used. This version implements both versions of short-range damping that appear in the literature (namely zero-damping and Becke-Johnson damping) provided the required parameters are specified manually, or the density functional is recognized from the input file - in which case default values will be used. This program was developed to analyze interatomic and intermolecular dispersion energies within the D3-scheme: if two molecules are recognized based on the interatomic connectivity then it is possible to ignore intramolecular terms.

If a density functional is not recognizable from the input/output file it will be necessary to specify the desired damping parameters. For zero-damping three terms (S6, S8 and RS6) are required. For Becke-Johnson damping four (S6, S8, A1, A2) are required.

##### Usage
```python
python dftd3.py (-damp zero/bj) (-s6 val) (-rs6 val) (-s8 val) (-a1 val) (-a2 val) (-im on/off) (-pw on/off) file(s)
```

The 3-body Axilrod-Teller-Muto 3-body dispersion terms can be switched on (they are not computed by default) with the flag ```-3body on```  

The type of damping can be switched from Becke-Johnson (the default) to Zero with the flag ```-damp zero```

To view the pairwise contribution of all dispersion terms ```-pw on``` (default is off)

To view only the intermolecular dispersion energy ```-im on``` (default is off, and requires the interatomic connectivity information to be specified in an input file so that there are two separate molecules)

#### Example 1:
Calculating the zero-damped D3 dispersion energy for an input file, manually specifying the relevant parameters for B3LYP.

```python
python dftd3.py 1.0 1.261 1.703 CH3F.com

s6 = 1.0 rs6 =  1.261 s8 = 1.703
   Breakdown   Attractive-R6   Attractive-R8   Repulsive-3-Body   Total   (Hartree)
   CH3F.com -1.92973198055e-05 -0.000722345327786 2.55309958043e-07 -0.000741387337633
```

#### Example 2:
Calculating the BJ-damped D3 dispersion energy for an output file, for which the functional (B3LYP) is detected and the parameters defined automatically.

```python
python dftd3.py CH3F2TS.log

D3-dispersion correction with Becke_Johnson damping: detected B3LYP functional - using default BJ-damping parameters
BJ-damping parameters: s6 = 1 s8 = 0.3981 a1 = 1.9889 a2 = 4.4211
3-body term will not be calculated

                                                 D3(R6)                   D3(R8)               Total (au)
examples/CH3F2TS.log              -0.000023517367523041    -0.000000976022647066    -0.000024493390170107
```

---
License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
