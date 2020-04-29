![pyDFTD3 logo](pyDFTD3_banner.png)

[![DOI](https://zenodo.org/badge/54939983.svg)](https://zenodo.org/badge/latestdoi/54939983)

This program will compute the Grimme D3-dispersion energy for a set of atomic Cartesian coordinates. This version implements both versions of short-range damping that appear in the literature (i.e., zero-damping and Becke-Johnson damping) provided the required parameters are specified manually, or the density functional can be automatically recognized from a Gaussian formatted output file - in which case default values will be used. This program was developed to analyze interatomic and intermolecular dispersion energies within the D3-scheme: if two molecules are recognized based on the interatomic connectivity then it is possible to ignore intramolecular terms.

If a density functional is not recognizable from the input/output file it will be necessary to specify the desired damping parameters. For zero-damping three terms (S6, S8 and RS6) are required. For Becke-Johnson damping four (S6, S8, A1, A2) are required.

The 3-body Axilrod-Teller-Muto 3-body dispersion terms can be switched on (they are not computed by default)

This program is no longer actively developed. Usage is now recommended through the [GoodVibes](https://github.com/bobbypaton/GoodVibes) program, where download of these files is necessary for this to function.

## Examples

Structure files available in the examples directory

1. D3-energy correction with zero-damping for a Gaussian output file. The density functional is parsed from the output and the appropriate damping parameters (s6, rs6, rs8) are applied automatically.

 ```python -m dftd3 formic_acid_dimer.log 
 
    D3-dispersion correction with zero-damping: detected B3LYP functional - using default zero-damping parameters
    Zero-damping parameters: s6 = 1.0 rs6 = 1.261 s8 = 1.703
    3-body term will not be calculated

                                      D3(R6)         D3(R8)         Total (au)
   formic_acid_dimer.log             -0.00088974    -0.50891629    -0.50980603
```

2. D3-energy correction with BJ-damping for a Gaussian output file. The density functional is parsed from the output and the appropriate damping parameters (s6, s8, a1, a2) are applied automatically.

 ```python -m dftd3 formic_acid_dimer.log -damp bj

    D3-dispersion correction with Becke_Johnson damping: detected B3LYP functional - using default BJ-damping parameters
    BJ-damping parameters: s6 = 1 s8 = 1.9889 a1 = 0.3981 a2 = 4.4211
    3-body term will not be calculated

                                      D3(R6)         D3(R8)         Total (au)
   formic_acid_dimer.log             -0.00455241    -0.00457708    -0.00912948
```

3. Pairwise breakdown of the D3(BJ) correction by atom pair.

 ```python -m dftd3 formic_acid_dimer.log -damp bj -pw

   --- Pairwise interaction between atoms 1 and 2 : Edisp = -0.359915 kcal/mol -0.3599152449643698
   --- Pairwise interaction between atoms 1 and 3 : Edisp = -0.320618 kcal/mol -0.6805336005703304
   ...
   --- Pairwise interaction between atoms 8 and 10 : Edisp = -0.136812 kcal/mol -5.6753050880809095
   --- Pairwise interaction between atoms 9 and 10 : Edisp = -0.053534 kcal/mol -5.728838870721404
```

4. D3-energy correction with BJ-damping for a Gaussian input file. The density functional is parsed from the route line. Based on the connectivity, we request only the intermolecular contributions to the D3-energy.

 ```python -m dftd3 formic_acid_dimer.com -damp bj -im

   Only computing intermolecular dispersion interactions! This is not the total D3-correction

   --- Ignoring interaction between atoms 1 and 2
   --- Ignoring interaction between atoms 1 and 3
   ...
   --- Ignoring interaction between atoms 8 and 10
   --- Ignoring interaction between atoms 9 and 10

                                      D3(R6)         D3(R8)         Total (au)
   formic_acid_dimer.com             -0.00152126    -0.00144858    -0.00296984
```

5. D3-energy correction with BJ-damping for an XYZ or PDB input file. Either the density functional or the parameters have to be specified in this case.

 ```python -m dftd3 formic_acid_dimer.xyz -damp bj -func b3lyp

   D3-dispersion correction with Becke_Johnson damping: detected B3LYP functional - using default BJ-damping parameters
   BJ-damping parameters: s6 = 1 s8 = 1.9889 a1 = 0.3981 a2 = 4.4211
   3-body term will not be calculated

                                      D3(R6)         D3(R8)         Total (au)
   formic_acid_dimer.xyz             -0.00455241    -0.00457708    -0.00912949

```

---
License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
