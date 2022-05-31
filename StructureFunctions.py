from EWVertex import *
from EWSE import *
from Charged_F import *


############ Parton structure functions ############

######## F1 - unpolarized ########
def F1A(Q2, x, MuR2, NLO):
    qi=1
    f1aph = 0
    if NLO == 0:
        while qi <= Nf:
            f1aph = f1aph + 0.25*hVq(Q2,MuR2,qi,0)**2*(PDFunpol(qi, x, MuR2)+PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f1aph
    if NLO == 1:
        while qi <= Nf:
            f1aph = (f1aph + (hAq(Q2,MuR2,qi,0)*hAq(Q2,MuR2,qi,1)
                    + hVq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1) + hVq(Q2,MuR2,qi,0)**2*SEAA(-Q2,MuR2)/Q2
                    -hVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*2*SW*CW*EtaAZ(Q2,1)*SEAZ(-Q2,MuR2)/Q2)*(PDFunpol(qi, x, MuR2)+PDFunpol(-qi, x, MuR2)))
            qi = qi + 1
        return f1aph

def F1AZ(Q2, x, MuR2, NLO):
    qi=1
    f1aph = 0
    if NLO == 0:
        while qi <= Nf:
            f1aph = f1aph + 0.5*(gVq(Q2,MuR2,qi,0)*hVq(Q2, MuR2, qi,0))*(PDFunpol(qi, x, MuR2)+PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f1aph

    if NLO == 1:
        while qi <= Nf:
            f1aph = (f1aph + ((gAq(Q2,MuR2,qi,1)*hAq(Q2,MuR2,qi,0) + gAq(Q2,MuR2,qi,0)*hAq(Q2,MuR2,qi,1)
                    + gVq(Q2,MuR2,qi,1)*hVq(Q2,MuR2,qi,0) + gVq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1))
                    + hVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*(SEAA(-Q2,MuR2)/Q2+SEZZ(-Q2,MuR2)/(Q2+MZ2))
                    - (hVq(Q2,MuR2,qi,0)**2 +EtaAZ(Q2,1)*(gVq(Q2,MuR2,qi,0)**2
                    + gAq(Q2,MuR2,qi,0)**2))*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFunpol(qi, x, MuR2)+PDFunpol(-qi, x, MuR2)))
            qi = qi + 1
        return f1aph

def F1Z(Q2, x, MuR2, NLO):
    qi=1
    f1aph = 0
    if NLO ==0:
        while qi <= Nf:
            f1aph = f1aph + 0.25*(gVq(Q2, MuR2, qi,0)**2 + gAq(Q2,qi, MuR2, 0)**2)*(PDFunpol(qi, x, MuR2)+PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f1aph

    if NLO == 1:
        while qi <= Nf:
            f1aph = (f1aph + (gAq(Q2, MuR2, qi,0)*gAq(Q2, MuR2, qi,1)
                    + gVq(Q2, MuR2, qi,0)*gVq(Q2, MuR2, qi,1)
                    + (gVq(Q2, MuR2, qi,0)**2 + gAq(Q2, MuR2, qi,0)**2)*SEZZ(-Q2,MuR2)/(Q2+MZ2)
                    - hVq(Q2, MuR2, qi,0)*gVq(Q2, MuR2, qi,0)*EtaAZ(Q2,1)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFunpol(qi, x, MuR2)+PDFunpol(-qi, x, MuR2)))
            qi = qi + 1
        return f1aph



######## F3 - unpolarized ########

def F3A(Q2, x, MuR2, NLO):
    qi=1
    f3aph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        while qi <= Nf:
            f3aph = f3aph + (2*(hAq(Q2, MuR2, qi,1)*hVq(Q2, MuR2, qi,0)
                    + hAq(Q2, MuR2, qi,0)*hVq(Q2, MuR2, qi,1))
                    +2*hVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*2*SW*CW*EtaAZ(Q2,1)*SEAZ(-Q2,MuR2)/Q2)*(PDFunpol(qi, x, MuR2)-PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f3aph


def F3AZ(Q2, x, MuR2, NLO):
    qi=1
    f3aph = 0
    if NLO ==0:
        while qi <= Nf:
            f3aph = f3aph + hVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*(PDFunpol(qi, x, MuR2)-PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f3aph
    if NLO == 1:
        while qi <= Nf:
            f3aph = f3aph + (2*(gVq(Q2,MuR2,qi,1)*hAq(Q2,MuR2,qi,0) + gVq(Q2,MuR2,qi,0)*hAq(Q2,MuR2,qi,1)
                    + gAq(Q2,MuR2,qi,1)*hVq(Q2,MuR2,qi,0)
                    + gAq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1))
                    + 2*hVq(Q2, MuR2, qi,0)*gAq(Q2, MuR2, qi,0)*(SEAA(-Q2,MuR2)/Q2+SEZZ(Q2,MuR2)/(Q2+MZ2))
                    - gAq(Q2, MuR2, qi,0)*gVq(Q2, MuR2, qi,0)*EtaAZ(Q2,1)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFunpol(qi, x, MuR2)-PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f3aph

def F3Z(Q2, x, MuR2, NLO):
    qi=1
    f3aph = 0
    if NLO ==0:
        while qi <= Nf:
            f3aph = f3aph + gVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*(PDFunpol(qi, x, MuR2)-PDFunpol(-qi, x, MuR2))
            qi = qi + 1
        return f3aph
    if NLO == 1:
        while qi <= Nf:
            f3aph = (f3aph +(2*(gAq(Q2,MuR2,qi,1)*gVq(Q2,MuR2,qi,0)
                    + gAq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,1)) + 4*gVq(Q2, MuR2, qi,0)*gAq(Q2, MuR2, qi,0)*SEZZ(-Q2,MuR2)/(Q2+MZ2)
                    - 2*hVq(Q2, MuR2, qi,0)*gAq(Q2, MuR2, qi,0)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFunpol(qi, x, MuR2)-PDFunpol(-qi, x, MuR2)))
            qi = qi + 1
        return f3aph


######## FL - unpolarized ########

def FLA(Q2, x, MuR2, NLO):
    qi=1
    flaph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        return 0

def FLAZ(Q2, x, MuR2, NLO):
    qi=1
    flaph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        return 0


def FLZ(Q2, x, MuR2, NLO):
    qi=1
    flaph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        return 0

######## g1 - polarized ########

def g1A(Q2, x, MuR2, NLO):
    qi=1
    g1aph = 0
    if NLO ==0:
        while qi <= Nf:
            g1aph = g1aph + hVq(Q2,MuR2,qi,0)**2/2*(PDFpol(qi, x, MuR2)+PDFpol(-qi, x, MuR2))
            qi = qi + 1
        return g1aph
    if NLO == 1:
        while qi <= Nf:
            g1aph = (g1aph + (hAq(Q2,MuR2,qi,0)*hAq(Q2,MuR2,qi,1)
                    + hVq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1)+hVq(Q2,MuR2,qi,0)**2*SEAA(-Q2,MuR2)/Q2
                    - hVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*2*SW*CW*EtaAZ(Q2,1)*SEAZ(-Q2,MuR2)/Q2) *(PDFpol(qi, x, MuR2)+PDFpol(-qi, x, MuR2)))
            qi = qi + 1
        return g1aph


def g1AZ(Q2, x, MuR2, NLO):
    qi=1
    g1aph = 0
    if NLO ==0:
        while qi <= Nf:
            g1aph = g1aph + hVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*(PDFpol(qi, x, MuR2)+PDFpol(-qi, x, MuR2))
            qi = qi + 1
        return g1aph
    if NLO == 1:
        while qi <= Nf:
            g1aph = (g1aph + (gAq(Q2,MuR2,qi,1)*hAq(Q2,MuR2,qi,0) + gAq(Q2,MuR2,qi,0)*hAq(Q2,MuR2,qi,1)
                             + gVq(Q2,MuR2,qi,1)*hVq(Q2,MuR2,qi,0)
                             + gVq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1)
                            + hVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*(SEAA(-Q2,MuR2)/Q2+SEZZ(-Q2,MuR2)/(Q2+MZ2))
                            - (hVq(Q2,MuR2,qi,0)**2 + EtaAZ(Q2,1)*(gAq(Q2,MuR2,qi,0)**2
                            +gVq(Q2,MuR2,qi,0)**2)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2))*(PDFpol(qi, x, MuR2)+PDFpol(-qi, x, MuR2)))
            qi = qi + 1
        return g1aph

def g1Z(Q2, x, MuR2, NLO):
    qi=1
    g1aph = 0
    if NLO ==0:
        while qi <= Nf:
            g1aph = g1aph + 0.5*(gVq(Q2,MuR2,qi,0)**2 + gAq(Q2,MuR2,qi,0)**2)*(PDFpol(qi, x, MuR2)+PDFpol(-qi, x, MuR2))
            qi = qi + 1
        return g1aph
    if NLO == 1:
        while qi <= Nf:
            g1aph = (g1aph + (gAq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,1)
                    + gVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,1) + (gAq(Q2,MuR2,qi,0)**2 + gVq(Q2,MuR2,qi,0)**2)*SEZZ(-Q2,MuR2)/(Q2+MZ2)
                    - hVq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFpol(qi, x, MuR2)+PDFpol(-qi, x, MuR2)))
            qi = qi + 1
        return g1aph


######## g5 - polarized ########

def g5A(Q2, x, MuR2, NLO):
    qi=1
    g5aph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        while qi <= Nf:
            g5aph = (g5aph +(hAq(Q2,MuR2,qi,1)*hVq(Q2,MuR2,qi,0)
                    + hAq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1)
                    - hVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*2*SW*CW*EtaAZ(Q2,1)*SEAZ(-Q2,MuR2)/Q2)*(PDFpol(qi, x, MuR2)-PDFpol(-qi, x, MuR2)))
            qi = qi + 1
        return g5aph


def g5AZ(Q2, x, MuR2, NLO):
    qi=1
    g5aph = 0
    if NLO == 0:
        while qi <= Nf:
            g5aph = g5aph - hVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*(PDFpol(qi, x, MuR2)-PDFpol(-qi, x, MuR2))
            qi = qi + 1
        return g5aph
    if NLO == 1:
        while qi <= Nf:
            g5aph = (g5aph + (gVq(Q2,MuR2,qi,1)*hAq(Q2,MuR2,qi,0) + gVq(Q2,MuR2,qi,0)*hAq(Q2,MuR2,qi,1)
                    + gAq(Q2,MuR2,qi,1)*hVq(Q2,MuR2,qi,0)
                    + gAq(Q2,MuR2,qi,0)*hVq(Q2,MuR2,qi,1)
                    + hVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*(SEAA(-Q2,MuR2)/Q2+SEZZ(-Q2,MuR2)/(Q2+MZ2))
                    - 2*gAq(Q2,MuR2,qi,0)*gVq(Q2,MuR2,qi,0)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFpol(qi, x, MuR2)-PDFpol(-qi, x, MuR2)))
            qi = qi + 1
        return g5aph

def g5Z(Q2, x, MuR2, NLO):
    qi=1
    g5aph = 0
    if NLO ==0:
        while qi <= Nf:
            g5aph = g5aph - gVq(Q2,MuR2, qi,0)*gAq(Q2,MuR2,qi,0)*(PDFpol(qi, x, MuR2)-PDFpol(-qi, x, MuR2))
            qi = qi + 1
        return g5aph
    if NLO == 1:
        while qi <= Nf:
            g5aph = (g5aph + (gAq(Q2,MuR2, qi,1)*gVq(Q2,MuR2, qi,0)
                    + gAq(Q2,MuR2, qi,0)*gVq(Q2,MuR2, qi,1)+2*gVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*SEZZ(-Q2,MuR2)/(Q2+MZ2)
                    -hVq(Q2,MuR2,qi,0)*gAq(Q2,MuR2,qi,0)*2*SW*CW*SEAZ(-Q2,MuR2)/Q2)*(PDFpol(qi, x, MuR2)-PDFpol(-qi, x, MuR2)))
            qi = qi + 1
        return g5aph


######## gL - polarized ########

def gLA(Q2, x, MuR2, NLO):
    qi=1
    glaph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        return 0

def gLAZ(Q2, x, MuR2, NLO):
    qi=1
    glaph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        return 0

def gLZ(Q2, x, MuR2, NLO):
    qi=1
    glaph = 0
    if NLO ==0:
        return 0
    if NLO == 1:
        return 0

################ Total structure functions assembled ################

######## Unpolarized Functions ########

def F1tot(Q2, x, MuR2, NLO):

    # charged current
    if currid ==0:
        qi = 1
        F1test = 0
        if NLO == 0:
            while qi <= (Nf+1)/2:
                F1test = F1test + PDFunpol(2*qi*LepType, x, MuR2) + PDFunpol(-(2*qi-1)*LepType, x, MuR2)
                qi = qi+1
            return F1test
        if NLO == 1:
            while qi <= (Nf+1)/2:
                F1test = F1test + (PDFunpol(2*qi*LepType, x, MuR2) + PDFunpol(-(2*qi-1)*LepType, x, MuR2))*(F1euQuark(-Q2, MuR2)+F1euLEP(-Q2, MuR2))
                qi = qi+1
            return F1test



    # neutral current
    if currid ==1:
        if NLO == 0:
            F1test = ((EtaAZ(Q2,2)*F1Z(Q2,x,MuR2,0)*(gAe(Q2, MuR2, 0)**2 + gVe(Q2, MuR2, 0)**2 - 2*gAe(Q2, MuR2, 0)*gVe(Q2, MuR2, 0)*LamE) +
                      F1A(Q2,x,MuR2,0)*(hAe(Q2, MuR2, 0)**2 + hVe(Q2, MuR2, 0)**2 - 2*hAe(Q2, MuR2, 0)*hVe(Q2, MuR2, 0)*LamE) +
                      EtaAZ(Q2,1)*F1AZ(Q2,x,MuR2,0)*(gVe(Q2, MuR2, 0)*(hVe(Q2, MuR2, 0) - hAe(Q2, MuR2, 0)*LamE) + gAe(Q2, MuR2, 0)*(hAe(Q2, MuR2, 0)- hVe(Q2, MuR2, 0)*LamE))))
            return F1test
        if NLO == 1:
            F1test = (F1A(Q2, x, MuR2, 1)*(hAe(Q2, MuR2, 0)**2 + hVe(Q2, MuR2, 0)**2 - 2*hAe(Q2, MuR2, 0)*hVe(Q2, MuR2, 0)*LamE) +
                      F1A(Q2, x, MuR2, 0)*(2*hAe(Q2, MuR2, 0)*hAe(Q2, MuR2, 1) + 2*hVe(Q2, MuR2, 0)*hVe(Q2, MuR2, 1) - 2*(hAe(Q2, MuR2, 1)*hVe(Q2, MuR2, 0) + hAe(Q2, MuR2, 0)*hVe(Q2, MuR2, 1))*LamE) +
                      EtaAZ(Q2, 1)*(F1Z(Q2, x, MuR2, 1)*(gAe(Q2, MuR2, 0)**2 + gVe(Q2, MuR2, 0)**2 - 2*gAe(Q2, MuR2, 0)*gVe(Q2, MuR2, 0)*LamE) +
                      F1Z(Q2,x,MuR2,0)*(2*gAe(Q2, MuR2, 0)*gAe(Q2, MuR2, 1) + 2*gVe(Q2, MuR2, 0)*gVe(Q2, MuR2, 1) -2*(gAe(Q2, MuR2, 1)*gVe(Q2, MuR2, 0) + gAe(Q2, MuR2, 0)*gVe(Q2, MuR2, 1))*LamE)) +
                      EtaAZ(Q2, 1)*(F1AZ(Q2,x,MuR2, 1)*(gVe(Q2, MuR2, 0)*(hVe(Q2, MuR2, 0) - hAe(Q2, MuR2, 0)*LamE) + gAe(Q2, MuR2, 0)*(hAe(Q2, MuR2, 0)- hVe(Q2, MuR2, 0)*LamE)) +
                      F1AZ(Q2, x, MuR2, 0)*(gVe(Q2, MuR2, 1)*(hVe(Q2, MuR2, 0) - hAe(Q2, MuR2, 0)*LamE) + gVe(Q2, MuR2, 0)*(hVe(Q2, MuR2, 1) - hAe(Q2, MuR2, 1)*LamE) +
                      gAe(Q2, MuR2, 1)*((hAe(Q2, MuR2, 0) - hVe(Q2, MuR2, 0)*LamE) + gAe(Q2, MuR2, 0)*(hAe(Q2, MuR2, 1) - hVe(Q2, MuR2, 1)*LamE)))))
            return F1test

def F3tot(Q2, x, MuR2, NLO):
    qi = 1
    F3test = 0

    # charged current
    if currid ==0:
        if NLO ==0:
            while qi <= (Nf+1)/2:
                F3test = F3test + 2*(PDFunpol(2*qi*LepType, x, MuR2) - PDFunpol(-(2*qi-1)*LepType, x, MuR2))*LepType
                qi = qi+1
            return F3test
        if NLO==1:
            while qi <= (Nf+1)/2:
                F3test = F3test + (PDFunpol(2*qi*LepType, x, MuR2) - PDFunpol(-(2*qi-1)*LepType, x, MuR2))*(F3euQuark(qi, -Q2, MuR2)+F3euLEP(qi,-Q2, MuR2))
                qi = qi+1
            return F3test

    # Neutral current
    if currid ==1:
        if NLO == 0:
            F3test = ((EtaAZ(Q2,2)*F3Z(Q2, x, MuR2, 0)*(2*gAe(Q2,MuR2,0)*gVe(Q2,MuR2,0) - (gAe(Q2,MuR2,0)**2 + gVe(Q2,MuR2,0)**2)*LamE) +
                    F3A(Q2,x,MuR2,0)*(2*hAe(Q2, MuR2, 0)*hVe(Q2,MuR2,0) - (hAe(Q2,MuR2,0)**2 + hVe(Q2,MuR2,0)**2)*LamE) +
                    EtaAZ(Q2,1)*F3AZ(Q2,x,MuR2,0)*(gAe(Q2,MuR2,0)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gVe(Q2,MuR2,0)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE))))
            return(F3test)
        if NLO == 1:
            F3test = (F3A(Q2,x,MuR2, 1)*(2*hAe(Q2,MuR2,0)*hVe(Q2,MuR2,0) - (hAe(Q2,MuR2,0)**2 + hVe(Q2,MuR2,0)**2)*LamE) +
                    2*F3A(Q2,x,MuR2,0)*(hAe(Q2,MuR2,1)*hVe(Q2,MuR2,0) + hAe(Q2,MuR2,0)*hVe(Q2,MuR2,1) - hAe(Q2,MuR2,0)*hAe(Q2,MuR2,1)*LamE - hVe(Q2,MuR2,0)*hVe(Q2,MuR2,1)*LamE) +
                    EtaAZ(Q2,2)*(F3Z(Q2,x,MuR2, 1)*(2*gAe(Q2,MuR2,0)*gVe(Q2,MuR2,0) - (gAe(Q2,MuR2,0)**2 + gVe(Q2,MuR2,0)**2)*LamE) +
                    2*F3Z(Q2,x,MuR2,0)*(gAe(Q2,MuR2,1)*gVe(Q2,MuR2,0) + gAe(Q2,MuR2,0)*gVe(Q2,MuR2,1) - gAe(Q2,MuR2,0)*gAe(Q2,MuR2,1)*LamE - gVe(Q2,MuR2,0)*gVe(Q2,MuR2,1)*LamE)) +
                    EtaAZ(Q2,1)*(F3AZ(Q2,x,MuR2, 1)*(gAe(Q2,MuR2,0)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gVe(Q2,MuR2,0)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE)) +
                    F3AZ(Q2,x,MuR2,0)*(gAe(Q2,MuR2,1)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gAe(Q2,MuR2,0)*(hVe(Q2,MuR2,1) - hAe(Q2,MuR2,1)*LamE) + gVe(Q2,MuR2,1)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE) +
                    gVe(Q2,MuR2,0)*(hAe(Q2,MuR2,1) - hVe(Q2,MuR2,1)*LamE))))
            return F3test


def FLtot(Q2, x, MuR2, NLO):
    if NLO ==0:
        return 0
    if NLO ==1:
        return 0
    if NLO ==2:
        return 0

######## Polarized Functions ########

def g1tot(Q2, x, MuR2, NLO):

    # charged current
    if currid ==0:
        qi = 1
        g1test = 0
        while qi <= (Nf+1)/2:
            g1test = g1test + PDFpol(2*qi*LepType, x, MuR2) + PDFpol(-(2*qi-1)*LepType, x, MuR2)
            qi = qi+1
        return(g1test)

    #neutral current
    if currid ==1:
        if NLO == 0:
            g1test = (EtaAZ(Q2,2)*g1Z(Q2,x,MuR2,0)*(2*gAe(Q2,MuR2,0)*gVe(Q2,MuR2,0) - (gAe(Q2,MuR2,0)**2 + gVe(Q2,MuR2,0)**2)*LamE) +
                      g1A(Q2,x,MuR2,0)*(2*hAe(Q2,MuR2,0)*hVe(Q2,MuR2,0) - (hAe(Q2,MuR2,0)**2 + hVe(Q2,MuR2,0)**2)*LamE) +
                      EtaAZ(Q2,1)*g1AZ(Q2,x,MuR2,0)*(gAe(Q2,MuR2,0)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gVe(Q2,MuR2,0)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE)))
            return g1test
        if NLO == 1:
            g1test = (g1A(Q2,x,MuR2, 1)*(2*hAe(Q2,MuR2,0)*hVe(Q2,MuR2,0) - (hAe(Q2,MuR2,0)**2 + hVe(Q2,MuR2,0)**2)*LamE) +
                      2*g1A(Q2,x,MuR2,0)*(hAe(Q2,MuR2,1)*hVe(Q2,MuR2,0) + hAe(Q2,MuR2,0)*hVe(Q2,MuR2,1) - hAe(Q2,MuR2,0)*hAe(Q2,MuR2,1)*LamE - hVe(Q2,MuR2,0)*hVe(Q2,MuR2,1)*LamE) +
                      EtaAZ(Q2,2)*(g1Z(Q2,x,MuR2, 1)*(2*gAe(Q2,MuR2,0)*gVe(Q2,MuR2,0) - (gAe(Q2,MuR2,0)**2 + gVe(Q2,MuR2,0)**2)*LamE) +
                      2*g1Z(Q2,x,MuR2,0)*(gAe(Q2,MuR2,1)*gVe(Q2,MuR2,0) + gAe(Q2,MuR2,0)*gVe(Q2,MuR2,1) - gAe(Q2,MuR2,0)*gAe(Q2,MuR2,1)*LamE - gVe(Q2,MuR2,0)*gVe(Q2,MuR2,1)*LamE)) +
                      EtaAZ(Q2,1)*(g1AZ(Q2,x,MuR2, 1)*(gAe(Q2,MuR2,0)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gVe(Q2,MuR2,0)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE)) +
                      g1AZ(Q2,x,MuR2,0)*(gAe(Q2,MuR2,1)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gAe(Q2,MuR2,0)*(hVe(Q2,MuR2,1) - hAe(Q2,MuR2,1)*LamE) + gVe(Q2,MuR2,1)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE) +
                      gVe(Q2,MuR2,0)*(hAe(Q2,MuR2,1) - hVe(Q2,MuR2,1)*LamE))))
            return g1test



def g5tot(Q2, x, MuR2, NLO):
    qi = 1
    g5test = 0

    #charged current
    if currid ==0:
        while qi <= (Nf+1)/2:
            g5test = g5test + (PDFpol(-(2*qi-1)*LepType, x, MuR2) - PDFpol(2*qi*LepType, x, MuR2))*LepType
            qi = qi+1
        return(g5test)

    #neutral current
    if currid ==1:
        if NLO ==0:
            g5test = (EtaAZ(Q2,2)*g5Z(Q2,x,MuR2,0)*(gAe(Q2,MuR2,0)**2 + gVe(Q2,MuR2,0)**2 - 2*gAe(Q2,MuR2,0)*gVe(Q2,MuR2,0)*LamE) +
                      g5A(Q2,x,MuR2,0)*(hAe(Q2,MuR2,0)**2 + hVe(Q2,MuR2,0)**2 - 2*hAe(Q2,MuR2,0)*hVe(Q2,MuR2,0)*LamE) +
                      EtaAZ(Q2,1)*g5AZ(Q2,x,MuR2,0)*(gVe(Q2,MuR2,0)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gAe(Q2,MuR2,0)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE)))
            return g5test
        if NLO ==1:
            g5test = (g5A(Q2,x,MuR2, 1)*(hAe(Q2,MuR2,0)**2 + hVe(Q2,MuR2,0)**2 - 2*hAe(Q2,MuR2,0)*hVe(Q2,MuR2,0)*LamE) +
                      g5A(Q2,x,MuR2,0)*(2*hAe(Q2,MuR2,0)*hAe(Q2,MuR2,1) + 2*hVe(Q2,MuR2,0)*hVe(Q2,MuR2,1) - 2*(hAe(Q2,MuR2,1)*hVe(Q2,MuR2,0) + hAe(Q2,MuR2,0)*hVe(Q2,MuR2,1))*LamE) +
                      EtaAZ(Q2,2)*(g5Z(Q2,x,MuR2, 1)*(gAe(Q2,MuR2,0)**2 + gVe(Q2,MuR2,0)**2 - 2*gAe(Q2,MuR2,0)*gVe(Q2,MuR2,0)*LamE) +
                      g5Z(Q2,x,MuR2,0)*(2*gAe(Q2,MuR2,0)*gAe(Q2,MuR2,1) + 2*gVe(Q2,MuR2,0)*gVe(Q2,MuR2,1) - 2*(gAe(Q2,MuR2,1)*gVe(Q2,MuR2,0) + gAe(Q2,MuR2,0)*gVe(Q2,MuR2,1))*LamE)) +
                      EtaAZ(Q2,1)*(g5AZ(Q2, x, MuR2, 1)*(gVe(Q2,MuR2,0)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gAe(Q2,MuR2,0)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE)) +
                      g5AZ(Q2,x,MuR2,0)*(gVe(Q2, MuR2, 1)*(hVe(Q2,MuR2,0) - hAe(Q2,MuR2,0)*LamE) + gVe(Q2,MuR2,0)*(hVe(Q2,MuR2,1) - hAe(Q2,MuR2,1)*LamE) + gAe(Q2,MuR2,1)*(hAe(Q2,MuR2,0) - hVe(Q2,MuR2,0)*LamE) +
                      gAe(Q2,MuR2,0)*(hAe(Q2, MuR2, 1) - hVe(Q2,MuR2,1)*LamE))))
            return g5test

def gLtot(Q2, x, MuR2, NLO):
    if currid ==0:
        return 0
    if currid ==1:
        if NLO ==0:
            gLtest = 0
            return gLtest
        if NLO ==1:
            return 0
        if NLO ==2:
            gltest = 0
            return gltest
