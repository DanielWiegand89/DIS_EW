###############################################
# Here's where we change all input parameters #
#                                             #
###############################################

import math as math

#### electron (1) positron (-1) ####
LepType = 1
ei = 1

#### Proton (1) ####
HadType = 1

####  Degree of polarization of lepton and hadron
LamE = 0.7
LamH = 0.7

#### Neutral (1) or Charged (0) Current ####
currid = 1

#### Numbers of active flavors ####
Nf = 5

#### Center of mass energy ####
Elep = 27.5
Ehad = 920
Shat = 4*Elep*Ehad

#### for differential cross section Q2 and x value ####
Q2test = 50**2
xtest = 0.1

#### Standard Model Parameters ####
aem = 1/137.036
alphasMZ = 0.1181
MZ = 91.1876
MW = 80.379
MW2 = MW**2
MZ2 = MZ**2
MT2 = 172.76**2
MB2 = 4.65**2
MH2 = 125.1**2
DeltaAlpha = 0.05954
SW = math.sqrt(1-MW**2/MZ**2)
SW2 = 1-MW**2/MZ**2
CW = MW/MZ
CW2 = CW**2
NAZ = Q2test/(Q2test+MZ**2)*(1/(4*SW**2*CW**2))
Mu2 = 1.0
MuF2 = Mu2

CF = 4/3
Tf = 1/2
Pi = math.pi
CA = 3
GeVtoPb = 0.389379*10**9

Zeta2 = Pi**2/6
Zeta3 = 1.2020569



