from Counterterms import *

#### Renormalized Self Energies

def SEAA(p12,MuR2):
    SEAA1 = (8*MB2 + 32*MT2 - 36*MW2 - 32*p12 - 8*A0fin(MB2,MuR2) - 32*A0fin(MT2,MuR2) + 36*A0fin(MW2,MuR2) +
             76*p12*B0fin(p12, 0, 0,MuR2) + 8*MB2*B0fin(p12, MB2, MB2,MuR2) + 4*p12*B0fin(p12, MB2, MB2,MuR2) +
             32*MT2*B0fin(p12, MT2, MT2,MuR2) + 16*p12*B0fin(p12, MT2, MT2,MuR2) -
             36*MW2*B0fin(p12, MW2, MW2,MuR2) - 27*p12*B0fin(p12, MW2, MW2,MuR2))
    return aem/(36*Pi)*SEAA1 + dZAA1(MuR2)*p12

def SEAZ(p12,MuR2):
    SEAZ1 = (2*(30*CW2*MW2 + 12*p12 + CW2*p12 - 6*MW2*SW2 - 31*p12*SW2 + MB2*(-6 + 8*SW2) +
            4*MT2*(-3 + 8*SW2)) - 4*(-3 + 4*SW2)*A0fin(MB2,MuR2) - 8*(-3 + 8*SW2)*A0fin(MT2,MuR2) +
             12*(-5*CW2 + SW2)*A0fin(MW2,MuR2) + 2*p12*(-27 + 76*SW2)*B0fin(p12, 0, 0,MuR2) +
             2*(2*MB2 + p12)*(-3 + 4*SW2)*B0fin(p12, MB2, MB2,MuR2) +
             4*(2*MT2 + p12)*(-3 + 8*SW2)*B0fin(p12, MT2, MT2,MuR2) +
             3*(32*CW2*MW2 + 19*CW2*p12 + 8*MW2*SW2 + p12*SW2)*B0fin(p12, MW2, MW2,MuR2))
    return aem/(72*CW*Pi*SW)*SEAZ1 + dZZA1(MuR2)/2*(p12-MZ2) + p12/2*dZAZ1(MuR2)

def SEZZ(p12,MuR2):
    SEZZ1 = (-2*CW2*(3*MH2 - 18*MT2 + 54*CW2**2*MW2 + 3*MZ2 + 23*p12 + 3*CW2**2*p12 + 48*MT2*SW2 -
            12*CW2*MW2*SW2 - 48*p12*SW2 + 2*CW2*p12*SW2 - 64*MT2*SW2**2 + 6*MW2*SW2**2 +
            63*p12*SW2**2 - 2*MB2*(9 - 12*SW2 + 8*SW2**2)) - 4*CW2*(9 - 12*SW2 + 8*SW2**2)*
             A0fin(MB2,MuR2) + (3*CW2*(-MH2 + MZ2 + 2*p12)*A0fin(MH2,MuR2))/p12 -
             4*CW2*(9 - 24*SW2 + 32*SW2**2)*A0fin(MT2,MuR2) + 12*CW2*(9*CW2**2 - 2*CW2*SW2 + SW2**2)*
             A0fin(MW2,MuR2) + (3*CW2*(MH2 - MZ2 + 2*p12)*A0fin(MZ2,MuR2))/p12 +
             4*CW2*p12*(27 - 54*SW2 + 76*SW2**2)*B0fin(p12, 0, 0,MuR2) +
             2*CW2*(p12*(9 - 12*SW2 + 8*SW2**2) + MB2*(-9 - 24*SW2 + 16*SW2**2))*
             B0fin(p12, MB2, MB2,MuR2) + (3*(12*MW2*p12 + CW2*(MH2**2 +
             (MZ2 - p12)**2 - 2*MH2*(MZ2 + p12)))*B0fin(p12, MH2, MZ2,MuR2))/p12
             + 2*CW2*(p12*(9 - 24*SW2 + 32*SW2**2) + MT2*(-9 - 48*SW2 + 64*SW2**2))*
             B0fin(p12, MT2, MT2,MuR2) - 3*CW2*(CW2**2*(60*MW2 + 39*p12) + 2*CW2*(-4*MW2 + p12)*SW2 -
            (20*MW2 + p12)*SW2**2)*B0fin(p12, MW2, MW2,MuR2))
    return aem/(144*CW2**2*Pi*SW2)*SEZZ1 +dZZZ1(MuR2)*(p12-MZ2)-dMZsq1(MuR2)

def SEWW(p12,Mu2):
    SEWW1 = (-(CW2*(-18*MB2**2 + 3*MH2**2 + 36*MB2*MT2 - 18*MT2**2 - 6*MH2*MW2 + 6*MW2**2 +
            24*CW2*MW2**2 - 6*MW2*MZ2 - 48*CW2*MW2*MZ2 + 3*MZ2**2 + 24*CW2*MZ2**2 -
            36*MW2*p12 + 36*CW2*MW2*p12 + 44*p12**2 + 8*CW2*p12**2 + 24*MW2**2*SW2 +
            36*MW2*p12*SW2 + 8*p12**2*SW2)) + 108*CW2*p12**2*B0fin(p12, 0, 0, Mu2) +
             24*CW2*(MW2**2 - 2*MW2*p12 - 5*p12**2)*SW2*B0fin(p12, 0, MW2, Mu2) -
             18*CW2*(MB2**2 + MT2**2 + MT2*p12 - 2*p12**2 + MB2*(-2*MT2 + p12))*
             B0fin(p12, MB2, MT2, Mu2) + 3*CW2*(MH2**2 + MW2**2 + 10*MW2*p12 + p12**2 -
            2*MH2*(MW2 + p12))*B0fin(p12, MH2, MW2, Mu2) +
             3*(CW2*(MW2**2 + (MZ2 - p12)**2 - 2*MW2*(MZ2 + p12)) +
            4*CW2**2*(2*MW2**2 + 2*MZ2**2 - 7*MZ2*p12 - 10*p12**2 - MW2*(4*MZ2 + 7*p12)) +
            12*MW2*p12*SW2**2)*B0fin(p12, MW2, MZ2, Mu2) + 18*CW2*MB2*(-MB2 + MT2 + 2*p12)*
            cm.log(MB2/Mu2) + 3*CW2*MH2*(MH2 - MW2 - 2*p12)*cm.log(MH2/Mu2) +
             18*CW2*MT2*(MB2 - MT2 + 2*p12)*cm.log(MT2/Mu2) +
             3*CW2*MW2*(-MH2 - MZ2 - 8*CW2*MZ2 - 40*p12 + 20*CW2*p12 + 20*p12*SW2 +
            MW2*(2 + 8*CW2 + 8*SW2))*cm.log(MW2/Mu2) + 3*CW2*(1 + 8*CW2)*MZ2*
             (-MW2 + MZ2 - 2*p12)*cm.log(MZ2/Mu2))
    return aem/(16*9*CW2*p12*Pi*SW2)*SEWW1
