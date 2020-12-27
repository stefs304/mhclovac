"""
HOPT810101 - Hydrophilicity value (Hopp-Woods, 1981) - 0.27901503759398494
ZIMJ680104 - Isoelectric point (Zimmerman et al., 1968) - 0.26629824561403503
KARS160105 - Average eccentricity (Karkbara-Knisley, 2016) - 0.24662155388471177
TANS770107 - Normalized frequency of left-handed helix (Tanaka-Scheraga, 1977) - 0.23065664160401
ISOY800104 - Normalized relative frequency of bend R (Isogai et al., 1980) - 0.2224235588972431
VELV850101 - Electron-ion interaction potential (Veljkovic et al., 1985) - 0.2001127819548872
MAXF760103 - Normalized frequency of zeta R (Maxfield-Scheraga, 1976) - 0.19167418546365914
RACS820104 - Average relative fractional occurrence in EL(i) (Rackovsky-Scheraga, 1982) - 0.1827794486215539
CHAM830102 - A parameter defined from the residuals obtained from the best correlation of the Chou-Fasman parameter of beta-sheet (Charton-Charton, 1983) - 0.17092982456140351
"""


class Config:

    INDEX_ID_LIST = ['PRAM900101', 'FASG760101', 'ZIMJ680104', 'CHOP780201', 'PRAM820103', 'RACS820112', 'ROBB760107']
    SIGMA = 0.8
    OVERLAP_DISTANCE = 1
    N_DISCRETE_POINTS = 10

