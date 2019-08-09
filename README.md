![pyDFTD3 logo](pyDFTD3_banner.png)

This program will compute the Grimme D3-dispersion energy for a set of atomic Cartesian coordinates. This version implements both versions of short-range damping that appear in the literature (i.e., zero-damping and Becke-Johnson damping) provided the required parameters are specified manually, or the density functional can be automatically recognized from a Gaussian formatted output file - in which case default values will be used. This program was developed to analyze interatomic and intermolecular dispersion energies within the D3-scheme: if two molecules are recognized based on the interatomic connectivity then it is possible to ignore intramolecular terms.

If a density functional is not recognizable from the input/output file it will be necessary to specify the desired damping parameters. For zero-damping three terms (S6, S8 and RS6) are required. For Becke-Johnson damping four (S6, S8, A1, A2) are required.

The 3-body Axilrod-Teller-Muto 3-body dispersion terms can be switched on (they are not computed by default)

This program is no longer actively developed. Usage is now recommended through the [GoodVibes](https://github.com/bobbypaton/GoodVibes) program, where download of these files is necessary for this to function.

---
License: [CC-BY](https://creativecommons.org/licenses/by/3.0/)
