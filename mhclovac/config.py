"""
ROSM880102 - Side chain hydropathy, corrected for solvation (Roseman, 1988) - 0.28197135416666663
ZIMJ680104 - Isoelectric point (Zimmerman et al., 1968) - 0.2740677083333333
OOBM770104 - Average non-bonded energy per residue (Oobatake-Ooi, 1977) - 0.24981510416666666
SNEP660101 - Principal component I (Sneath, 1966) - 0.24546093750000003
ROBB760111 - Information measure for C-terminal turn (Robson-Suzuki, 1976) - 0.23538281249999998
CHAM830102 - A parameter defined from the residuals obtained from the best correlation of the Chou-Fasman parameter of beta-sheet (Charton-Charton, 1983) - 0.20620052083333332
RACS820104 - Average relative fractional occurrence in EL(i) (Rackovsky-Scheraga, 1982) - 0.18330468749999998
WERD780103 - Free energy change of alpha(Ri) to alpha(Rh) (Wertz-Scheraga, 1978) - 0.17045833333333335
KARS160120 - Weighted minimum eigenvalue based on the atomic numbers (Karkbara-Knisley, 2016) - 0.1413723958333333

"""


class Config:

    INDEX_ID_LIST = ['ROSM880102', 'ZIMJ680104', 'OOBM770104', 'SNEP660101', 'ROBB760111', 'CHAM830102', 'RACS820104', 'WERD780103', 'KARS160120']
    SIGMA = 0.8
    OVERLAP_DISTANCE = 1
    N_DISCRETE_POINTS = 10

