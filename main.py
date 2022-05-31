###############################################
#            DEEP INELASTIC SCATTERING        #
#            NLO QCD, NLO EW, NNLO QCD        #
#                                             #
#  D. Wiegand                                 #
#  Contact: Daniel.Wiegand@northwestern.edu   #
#  Reference: arXiv                           #
###############################################

# QCD IS STILL WORK IN PROGRESS

from StructureFunctions import *
from BoxCont import *

# Unpol CS dxdy

def UnpolCS(x, Q2, corrid):
    y = Q2/(Shat*x)
    dsig = 2*x*(1+(1-y)**2)*F1tot(Q2, x, Mu2, corrid)
    dsig = dsig + LepType*x*(1-(1-y)**2)*F3tot(Q2, x, Mu2, corrid)
    dsig = dsig + 2*(1-y)*FLtot(Q2, x, Mu2, corrid)
    dsig = (dsig*Omegafct(Q2, currid)*2*Pi*aem**2/(x*y*Q2) + dsigmaBox(x, Q2, Mu2, 0, corrid))*GeVtoPb
    return dsig

# Pol CS dxdy

def PolCS(x, Q2, corrid):
    y = Q2/(Shat*x)
    dsig = x*(1+(1-y)**2)*g5tot(Q2, x, Mu2, corrid)
    dsig = dsig - LepType*x*(1-(1-y)**2)*g1tot(Q2, x, Mu2, corrid)
    dsig = dsig + (1-y)*gLtot(Q2, x, Mu2, corrid)
    dsig = (dsig*Omegafct(Q2, currid)*4*Pi*aem**2/(x*y*Q2) + dsigmaBox(x, Q2, Mu2, 1, corrid))*GeVtoPb
    return dsig

# CS with incomplete polarization

dsigfin0 = (1-LamH)*UnpolCS(xtest, Q2test, 0) + LamH*PolCS(xtest, Q2test, 0)

dsigfin1 = (1-LamH)*UnpolCS(xtest, Q2test, 1) + LamH*PolCS(xtest, Q2test, 1)

print("Cross Section", dsigfin0, "pb (LO), plus", np.real(dsigfin1), "pb")
