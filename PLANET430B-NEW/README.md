# Planet430B-New
It is based on highly precise SOFA algorithms for computing Solar System Bodies through the new model of precession-nutation adopted since 2009 Jan.01 by IAU (International Astronomy Union) in order to take into account all movements of Terrestrial Pole, detailed up to tidal and relativistic perturbations.

For such a reason the two databases shown in figure <b>updating_2dbases.png</b> should be periodically monitored and updated. Then a new compilation is to be performed and related (a.exe, for Windows systems) formed.

Find below the links to download both, coming from this URL: https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html
### 1) final-iers.txt 
Paste&copy from finals.data (IAU2000), latest version
https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt
Be careful: finals-iers.txt first row should start as follows:
92 1 1 48622.00 I  0.182987 0.000672 ....
### 2) iers1962_now.txt 
Paste&copy from EOP 14 C04 (IAU2000A) - one file (1962-now), latest version
https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt
Be careful: iers1962_now.txt first row should start as follows:
1962   1   1  37665  -0.012700   0.213000  ....

## IERS updating 
Fortran code lines 137,139 adjusted as follows

    FINAL_ = 59070       ! Aug.04,2019 (file final-iers.txt) 
    EOP_final = 58666    ! Aug.04,2019 (file iers1962_now.txt)
