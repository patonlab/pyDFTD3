%chk=CH3F2TS.chk
# freq rb3lyp/3-21g scrf=check guess=tcheck geom=(connectivity,allcheck)
genchk

Title Card Required

-1 1
 C                  0.00000000    0.00000000    0.00000000
 H                  0.00000000    1.07977500    0.00000000
 H                 -0.93511300   -0.53988800    0.00000000
 H                  0.93511300   -0.53988800    0.00000000
 F                  0.00000000    0.00000000    1.80755300
 F                  0.00000000    0.00000000   -1.80755300

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2
 3
 4
 5
 6

