from Input import *
from mpmath import polylog
import cmath as cm
import numpy as np
from scipy.integrate import quad

from AlphaS import *
from pdfflow import mkPDFs
import tensorflow as tf

# Load the PDFs we're using

pdfunpoltot = mkPDFs("NNPDF31_nlo_pdfas", [0], dirname="/Users/Mac/Downloads/ManeParse_2.0/Demo/PDF_Sets/LHA/")
pdfpoltot = mkPDFs("NNPDFpol11_100", [0], dirname="/Users/Mac/Downloads/ManeParse_2.0/Demo/PDF_Sets/LHA/")

def PDFunpol(pid, x1, Q2):
    if abs(pid) == 6:
        return 0
    else:
        pdftest = pdfunpoltot.py_xfxQ2([pid], [x1], [Q2])
        pdftest = pdftest.numpy()/x1
        return pdftest

def PDFpol(pid, x1, Q2):
    if abs(pid) == 6:
        return 0
    else:
        pdftest = pdfpoltot.py_xfxQ2([pid], [x1], [Q2])
        pdftest = pdftest.numpy()/x1
        return pdftest


################### Functions necessary to evaluate the loop integrals

def Omegafct(Q2, curr):
    if curr == 1:
        return 1.0
    if curr == 0:
        return (1-LamE*LepType)*(Q2/(Q2+MW2)*1/(64*MW2*SW2))**2

def EtaAZ(Q2, id):
    etaAZ = Q2/(Q2+MZ2)*(1/(4*SW2*CW2))
    if id ==1:
        return etaAZ
    if id ==2:
        return etaAZ**2

def li2(x):
    return polylog(2,x)

def li3(x):
    return polylog(3,x)

def S12(z):
    def S12int(t):
        return (np.real(cm.log(1-z*t))**2)/t
    S12res, S12err = quad(S12int, 0, 1)
    return S12res/2

################### Tadpole

def A0fin(m12,MuR2):
    return m12*(1-cm.log(m12/MuR2))

################### Bubble

def B0fin(p12,m12,m22,MuR2):
    if m12*m22 !=0:
        num1 = -m12 + m22 + p12 - cm.sqrt(m12**2 + (m22 - p12)**2 - 2*m12*(m22 + p12))
        num2 = -m12 + m22 + p12 + cm.sqrt(m12**2 + (m22 - p12)**2 - 2*m12*(m22 + p12))
        B0test = num1*cm.log(1 - 2*p12/num1) + num2*cm.log(1 - 2*p12/num2)
        B0test = B0test/(2*p12) + 2 - cm.log(m12/MuR2)
        return B0test
    if (m12==0 and m22==0):
        return 2-cm.log(p12/MuR2)
    else:
        M12 = m12+m22
        return 2 - (p12 - M12)/p12*cm.log(1 - p12/M12) - cm.log(M12/MuR2)

################### Triangles

def C01(S,m12):
    C0test = cm.log((m12 + S)/m12)*cm.log((m12 + S)/S)
    C0test = C0test + li2(m12/(m12 + S)) - li2((m12 + S)/m12) + li2(S/(m12 + S))
    return C0test/S

def C02(S,m12,m22):
    a1 =2*m12/(m12 + m22 - S - cm.sqrt(m12**2 + (m22 - S)**2 - 2*m12*(m22 + S)))
    a2 =2*(m12 - S)/(m12 + m22 - S - cm.sqrt(m12**2 + (m22 - S)**2 - 2*m12*(m22 + S)))
    a3 =2*m12/(m12 + m22 - S + cm.sqrt(m12**2 + (m22 - S)**2 - 2*m12*(m22 + S)))
    a4 =2*(m12 - S)/(m12 + m22 - S + cm.sqrt(m12**2 + (m22 - S)**2 - 2*m12*(m22 + S)))
    C0test = cm.log(S/m12)*cm.log(1 - S/m12) + li2(S/m12)
    C0test = C0test -li2(a1) +li2(a2) -li2(a3) + li2(a4)
    return -C0test/S

###### Boxes (TBD)

def D01(S,T,m12,m22):
    def intfct(x):
        gammas = cm.sqrt(m12**2 + (m22 + T*(-1 + x))**2 - 2*m12*(m22 + T - T*x))
        intfct = (2*cm.log(-(cm.sqrt(m12*m22)/(S*x))) - ((m12 + m22 - T + 2*S*x + T*x)*
                cm.log((gammas + m12 + m22 - T + T*x)/(-gammas + m12 + m22 - T + T*x)))/gammas)/(2*(m12*(m22 + S*x)
                                        + S*x*(m22 + T*(-1 + x) + S*x)))
        return np.real(intfct)

    D01res, D01err = quad(intfct, 0, 1)

    return D01res


################### Counter term functions ###################

# Boson mass renormalizations

def dMZsq1(MuR2):
    dMZsq1 = -aem/(144*CW2**2*MZ2*Pi*SW2)*(2*CW2*MZ2*(3*MH2 - 18*MT2 + 54*CW2**2*MW2 + 26*MZ2 + 3*CW2**2*MZ2 + 48*MT2*SW2 -
            12*CW2*MW2*SW2 - 48*MZ2*SW2 + 2*CW2*MZ2*SW2 - 64*MT2*SW2**2 + 6*MW2*SW2**2 +
            63*MZ2*SW2**2 - 2*MB2*(9 - 12*SW2 + 8*SW2**2)) + 4*CW2*MZ2*(9 - 12*SW2 + 8*SW2**2)*
            A0fin(MB2,MuR2) + 3*CW2*(MH2 - 3*MZ2)*A0fin(MH2,MuR2) + 4*CW2*MZ2*(9 - 24*SW2 + 32*SW2**2)*
            A0fin(MT2,MuR2) - 12*CW2*MZ2*(9*CW2**2 - 2*CW2*SW2 + SW2**2)*A0fin(MW2,MuR2) -
            3*CW2*(MH2 + MZ2)*A0fin(MZ2,MuR2) - 4*CW2*MZ2**2*(27 - 54*SW2 + 76*SW2**2)*B0fin(MZ2, 0, 0,MuR2) -
            2*CW2*MZ2*(MZ2*(9 - 12*SW2 + 8*SW2**2) + MB2*(-9 - 24*SW2 + 16*SW2**2))*
            B0fin(MZ2, MB2, MB2,MuR2) - 3*(CW2*MH2*(MH2 - 4*MZ2) + 12*MW2*MZ2)*B0fin(MZ2, MH2, MZ2,MuR2) -
            2*CW2*MZ2*(MZ2*(9 - 24*SW2 + 32*SW2**2) + MT2*(-9 - 48*SW2 + 64*SW2**2))*
            B0fin(MZ2, MT2, MT2,MuR2) - 3*CW2*MZ2*(-3*CW2**2*(20*MW2 + 13*MZ2) + (20*MW2 + MZ2)*SW2**2 +
            CW2*(8*MW2*SW2 - 2*MZ2*SW2))*B0fin(MZ2, MW2, MW2,MuR2))
    return dMZsq1

def dMWsq1(MuR2):
    dMWsq1 = -aem/(144*CW2*MW2*Pi*SW2)*(2*CW2*MW2*(-18*MB2 + 3*MH2 - 18*MT2 + 64*MW2 - 8*CW2*MW2 + 3*MZ2 + 24*CW2*MZ2 +
            64*MW2*SW2) - 18*CW2*(MB2 - MT2 - 2*MW2)*A0fin(MB2,MuR2) +
            3*CW2*(MH2 - 3*MW2)*A0fin(MH2,MuR2) + 18*CW2*(MB2 - MT2 + 2*MW2)*A0fin(MT2,MuR2) -
            3*CW2*(MH2 + MZ2 + 8*CW2*MZ2 + MW2*(38 - 28*CW2 - 76*SW2))*A0fin(MW2,MuR2) -
            3*CW2*(1 + 8*CW2)*(3*MW2 - MZ2)*A0fin(MZ2,MuR2) - 108*CW2*MW2**2*B0fin(MW2, 0, 0,MuR2) +
            18*CW2*((MB2 - MT2)**2 + (MB2 + MT2)*MW2 - 2*MW2**2)*B0fin(MW2, MB2, MT2,MuR2) -
            3*CW2*(MH2**2 - 4*MH2*MW2 + 12*MW2**2)*B0fin(MW2, MH2, MW2,MuR2) +
            (-3*CW2*MZ2*(-4*MW2 + MZ2) + 12*CW2**2*(15*MW2**2 + 11*MW2*MZ2 - 2*MZ2**2) -
            36*MW2**2*SW2**2)*B0fin(MW2, MW2, MZ2,MuR2))
    return dMWsq1

# Photon/Z field renormalization

def dZAA1(MuR2):
    dZAA1 = (-DeltaAlpha - (aem*(5*MT2*MW2 + 16*MW2*A0fin(MT2,MuR2) - 27*MT2*A0fin(MW2,MuR2)))/
        (36*MT2*MW2*Pi) - (aem*(6*MB2 - 20*MZ2 - 6*A0fin(MB2,MuR2) + 57*MZ2*B0fin(MZ2, 0, 0,MuR2) +
        6*MB2*B0fin(MZ2, MB2, MB2,MuR2) + 3*MZ2*B0fin(MZ2, MB2, MB2,MuR2)))/(27*MZ2*Pi))
    return dZAA1

def dZZA1(MuR2):
    dZZA1 = (aem*(-MW2 + A0fin(MW2,MuR2)))/(CW*MZ2*Pi*SW)
    return dZZA1

def dZAZ1(MuR2):
    dZAZ1 = aem/(36*CW*MZ2*Pi*SW)*(36*MW2 + 24*MW2*(3*CW2 + SW2) - 2*(CW2*(12*MW2 + 19*MZ2) + 12*MW2*SW2 +
            MZ2*(-24 + 65*SW2)) - 4*(-3 + 4*SW2)*(MB2 - A0fin(MB2,MuR2)) -
            8*(-3 + 8*SW2)*(MT2 - A0fin(MT2,MuR2)) - 12*(MW2*(9*CW2 - SW2) + (-5*CW2 + SW2)*A0fin(MW2,MuR2)) -
            2*MZ2*(-27 + 76*SW2)*(-1 + B0fin(MZ2, 0, 0,MuR2)) +
            (3 - 4*SW2)*(-2*MZ2 + 2*(2*MB2 + MZ2)*B0fin(MZ2, MB2, MB2,MuR2)) +
            2*(3 - 8*SW2)*(-2*MZ2 + 2*(2*MT2 + MZ2)*B0fin(MZ2, MT2, MT2,MuR2)) -
            3*(CW2*(32*MW2 + 19*MZ2) + (8*MW2 + MZ2)*SW2)*B0fin(MZ2, MW2, MW2,MuR2))
    return dZAZ1

def dZZZ1(MuR2):
    dZZZ1 = aem/(288*CW2**2*Pi*SW2)*((36*MW2)/(MH2 - MZ2) - (6*CW2**3*(120*MW2**2 + 10*MW2*MZ2 - 37*MZ2**2))/((4*MW2 - MZ2)*MZ2)
            + (4*CW2**2*(6*MW2 - MZ2)*SW2)/MZ2 +
            CW2*(236 - (12*MH2)/(MH2 - MZ2) + (252*MT2)/(4*MT2 - MZ2) +
            (3*MH2**2)/((MH2 - MZ2)*MZ2) - (72*MZ2)/(4*MB2 - MZ2) - (72*MZ2)/(4*MT2 - MZ2) -
            (72*MT2**2)/(4*MT2*MZ2 - MZ2**2) - 480*SW2 - (384*MT2*SW2)/(4*MT2 - MZ2) +
            (96*MZ2*SW2)/(4*MB2 - MZ2) + (192*MZ2*SW2)/(4*MT2 - MZ2) -
            (384*MT2**2*SW2)/(4*MT2*MZ2 - MZ2**2) + 700*SW2**2 + (512*MT2*SW2**2)/(4*MT2 - MZ2) -
            (108*MW2*SW2**2)/(4*MW2 - MZ2) - (64*MZ2*SW2**2)/(4*MB2 - MZ2) -
            (256*MZ2*SW2**2)/(4*MT2 - MZ2) - (6*MZ2*SW2**2)/(4*MW2 - MZ2) +
            (512*MT2**2*SW2**2)/(4*MT2*MZ2 - MZ2**2) + (240*MW2**2*SW2**2)/(4*MW2*MZ2 - MZ2**2) +
            (8*MB2**2*(-9 - 24*SW2 + 16*SW2**2))/((4*MB2 - MZ2)*MZ2) +
            (4*MB2*(63 - 48*SW2 + 32*SW2**2))/(4*MB2 - MZ2)) +
            (8*CW2*(MB2*(9 + 24*SW2 - 16*SW2**2) + MZ2*(-9 + 12*SW2 - 8*SW2**2))*A0fin(MB2,MuR2))/
            (4*MB2*MZ2 - MZ2**2) +
            (3*(-12*MW2*(MH2 - 2*MZ2)*MZ2 + CW2*(-3*MH2**3 + 12*MH2**2*MZ2 - 14*MH2*MZ2**2 + 2*MZ2**3))*
            A0fin(MH2,MuR2))/((MH2 - MZ2)**2*MZ2**2) +
            (8*CW2*(MT2*(9 + 48*SW2 - 64*SW2**2) + MZ2*(-9 + 24*SW2 - 32*SW2**2))*A0fin(MT2,MuR2))/
            (4*MT2*MZ2 - MZ2**2) + (12*CW2*(CW2**2*(60*MW2 + 39*MZ2) + 2*CW2*(-4*MW2 + MZ2)*SW2 -
            (20*MW2 + MZ2)*SW2**2)*A0fin(MW2,MuR2))/((4*MW2 - MZ2)*MZ2) +
            (3*(-12*MW2*MZ2**2 + CW2*(2*MH2**3 - 7*MH2**2*MZ2 + 10*MH2*MZ2**2 - 2*MZ2**3))*A0fin(MZ2,MuR2))/
            ((MH2 - MZ2)**2*MZ2**2) - 8*CW2*(27 - 54*SW2 + 76*SW2**2)*B0fin(MZ2, 0, 0,MuR2) +
            (4*CW2*(-2*MB2*MZ2*(9 - 12*SW2 + 8*SW2**2) + MZ2**2*(9 - 12*SW2 + 8*SW2**2) +
            2*MB2**2*(-9 - 24*SW2 + 16*SW2**2))*B0fin(MZ2, MB2, MB2,MuR2))/(4*MB2*MZ2 - MZ2**2) +
            (3*(CW2*MH2*(3*MH2 - 8*MZ2) + 12*MW2*MZ2)*B0fin(MZ2, MH2, MZ2,MuR2))/MZ2**2 +
            (4*CW2*(-2*MT2*MZ2*(9 - 24*SW2 + 32*SW2**2) + MZ2**2*(9 - 24*SW2 + 32*SW2**2) +
            2*MT2**2*(-9 - 48*SW2 + 64*SW2**2))*B0fin(MZ2, MT2, MT2,MuR2))/(4*MT2*MZ2 - MZ2**2) +
            (3*CW2*(6*CW2**2*(40*MW2**2 - 26*MW2*MZ2 + 13*MZ2**2) -
            4*CW2*(8*MW2**2 + 2*MW2*MZ2 - MZ2**2)*SW2 - 2*(40*MW2**2 - 2*MW2*MZ2 + MZ2**2)*SW2**2)*
            B0fin(MZ2, MW2, MW2,MuR2))/(-4*MW2*MZ2 + MZ2**2))
    return dZZZ1

def dZWW1(MuR2):
    dZWW1 = 1
    return dZWW1

# weak mixing angle renormalization

def dSW1(MuR2):
    dSW1 = CW2/(2*SW)*(dMZsq1(MuR2)/MZ2 - dMWsq1(MuR2)/MW2)
    return dSW1

# charge renormalization

def dZe1(MuR2):
    dZe1 = -0.5*dZAA1(MuR2) - SW/(2*CW)*dZZA1(MuR2)
    return dZe1

# external fermion wave function renormalizations

def dZnuL(MuR2):
    dZnuL = (aem*(1 + 2*CW2 + 4*CW2*cm.log(MW2/Mu2) + 2*cm.log(MZ2/Mu2)))/(32*CW2*Pi*SW2)
    return dZnuL

def dZeR(MuR2):
    dZfeR = (aem*SW2*(3*MZ2 - 2*A0fin(MZ2,MuR2)))/(8*MW2*Pi)
    return dZfeR

def dZeL(MuR2):
    dZfeL= (aem*(3*MW2*MZ2*(2*CW2 + (1 - 2*SW2)**2) - 4*CW2*MZ2*A0fin(MW2,MuR2)
            - 2*MW2*(1 - 2*SW2)**2*A0fin(MZ2,MuR2)))/(32*MW2**2*Pi*SW2)
    return dZfeL

def dZuR(MuR2):
    dZfuR = (aem*SW2*(3*MZ2 - 2*A0fin(MZ2,MuR2)))/(18*MW2*Pi)
    return dZfuR

def dZuL(MuR2):
    dZfuL = ((aem*(-36*CW2*MZ2*A0fin(MW2,MuR2) + MW2*(3*MZ2*(18*CW2 + (3 - 4*SW2)**2) -
            2*(3 - 4*SW2)**2*A0fin(MZ2,MuR2))))/(288*CW2*MW2*MZ2*Pi*SW2))
    return dZfuL

def dZdR(MuR2):
    dZfdR = (aem*SW2*(3*MZ2 - 2*A0fin(MZ2,MuR2)))/(72*MW2*Pi)
    return dZfdR

def dZdL(MuR2):
    dZfdL = (aem*(3*MW2*MZ2*(18*CW2 + (3 - 2*SW2)**2) - 36*CW2*MZ2*A0fin(MW2,MuR2)
            -2*MW2*(3 - 2*SW2)**2*A0fin(MZ2,MuR2)))/(288*MW2**2*Pi*SW2)
    return dZfdL



# W total counterterm

def dWct(Q2, MuR2):
    dWct = 4*dZe1(MuR2) - (2*dMWsq1(MuR2))/(MW2 + Q2) - (4*dSW1(MuR2))/SW + dZnuL(MuR2) +dZeL(MuR2) + dZuL(MuR2) + dZdL(MuR2)
    return dWct