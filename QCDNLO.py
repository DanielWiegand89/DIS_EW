# NLO QCD Corrections to the structure functions
# Checked against arXiv: 1210.7203

##### THIS IS STILL WORK IN PROGRESS

from AlphaS import *
from EWVertex import *
from pdfflow import mkPDFs


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

# NLO QCD unpolarized

def F1QCD(Q2, x, MuF2, q):

    def F1qintPP(z):
        F1qintPP = ((1/(1-z))*(2*cm.log(Q2/MuF2)-1.5) + 2*cm.log(1-z)/(1-z))
        F1qintPP = F1qintPP*((PDFunpol(q,x/z,MuF2)+PDFunpol(-q,x/z,MuF2))/z - PDFunpol(q,x,MuF2) - PDFunpol(-q,x,MuF2))
        return np.real(F1qintPP)

    F1qintPPres, F1qintPPerr = quad(F1qintPP, x, 1)

    def F1qintreg(z):
        F1qintreg = 1/z*(-(1+z)*cm.log(Q2/MuF2*(1-z))-(1+z**2)/(1-z)*cm.log(z)+3+2*z)*(PDFunpol(q,x/z,MuF2)+PDFunpol(-q,x/z,MuF2))
        return np.real(F1qintreg)

    F1qintregres, F1qintregerr = quad(F1qintreg, x, 1)


    F1qtot = F1qintPPres - (PDFunpol(q,x,Q2)+PDFunpol(-q,x,Q2))*(0.5*(3-4*cm.log(Q2/MuF2))*cm.log(1-x)-2*li2(x))
    F1qtot = CF*(F1qtot + F1qintregres + (1.5*cm.log(Q2/MuF2)-4.5-2*Zeta2)*(PDFunpol(q,x,Q2)-PDFunpol(-q,x,Q2)))

    def F1gint(z):
        F1gint = 1/z*0.5*((z**2+(1-z)**2)/2*(cm.log(Q2/MuF2*(1-z)/z)-1)+2*z*(1-z))*PDFunpol(0,x/z,MuF2)
        return np.real(F1gint)

    F1gintres, F1ginterr = quad(F1gint, x, 1)

    return (F1qtot + F1gintres)

def F3QCD(Q2, x, MuF2, q):

    def F3qintPP(z):
        F3qintPP = ((1/(1-z))*(2*cm.log(Q2/MuF2)-1.5) + 2*cm.log(1-z)/(1-z))
        F3qintPP = F3qintPP*((PDFunpol(q,x/z,MuF2)-PDFunpol(-q,x/z,MuF2))/z - PDFunpol(q,x,MuF2) + PDFunpol(-q,x,MuF2))
        return np.real(-F3qintPP)

    F3qintPPres, F3qintPPerr = quad(F3qintPP, x, 1)

    def F3qintreg(z):
        F3qintreg = (-(1+z)*cm.log(Q2/MuF2*(1-z))-(1+z**2)/(1-z)*cm.log(z)+2+3*z)*(PDFunpol(q,x/z,MuF2)-PDFunpol(-q,x/z,MuF2))
        return -np.real(F3qintreg)

    F3qintregres, F3qintregerr = quad(F3qintreg, x, 1)

    F3qtot = F3qintPPres - (PDFunpol(q,x,Q2)-PDFunpol(-q,x,Q2))*(0.5*(3-4*cm.log(Q2/MuF2))*cm.log(1-x)-2*li2(x))
    F3qtot = CF*(F3qtot + F3qintregres + (1.5*cm.log(Q2/MuF2)-4.5-2*Zeta2)*(PDFunpol(q,x,Q2)-PDFunpol(-q,x,Q2)))

    return als(MuF2, 2)/(2*Pi)*(F3qtot)


def FLQCD(Q2, x, MuF2, q):
    def FLqint(z):
        FLqint = 1/z*2*z*CF*(PDFunpol(q, x/z, MuF2)-PDFunpol(-q, x/z, MuF2))
        return np.real(FLqint)

    FLqintres, FLqinterr = quad(FLqint, x, 1)

    def FLgint(z):
        FLgint = 1/z*2*z*(1-z)*PDFunpol(0, x/z, MuF2)
        return np.real(FLgint)

    FLgintres, FLginterr = quad(FLgint, x, 1)

    FLAQCDtot = FLqintres + FLgintres

    return FLAQCDtot


#### NLO QCD Polarized

def g1QCD(Q2, x, MuF2, q):

    def g1qintPP(z):
        g1qintPP = ((1/(1-z))*(2*cm.log(Q2/MuF2)-1.5) + 2*cm.log(1-z)/(1-z))
        g1qintPP = g1qintPP*((PDFpol(q,x/z,MuF2)+PDFpol(-q,x/z,MuF2))/z - PDFpol(q,x,MuF2)- PDFpol(-q,x,MuF2))
        return np.real(g1qintPP)

    g1qintPPres, g1qintPPerr = quad(g1qintPP, x, 1)

    def g1qintreg(z):
        g1qintreg = (-(1+z)*cm.log(Q2/MuF2*(1-z))-(1+z**2)/(1-z)*cm.log(z)+2+3*z)*(PDFpol(q,x/z,MuF2)+PDFpol(-q,x/z,MuF2))
        return np.real(g1qintreg)

    g1qintregres, g1qintregerr = quad(g1qintreg, x, 1)

    g1qtot = g1qintPPres - (PDFpol(q,x,Q2)+PDFpol(-q,x,Q2))*(0.5*(3-4*cm.log(Q2/MuF2))*cm.log(1-x)-2*li2(x))
    g1qtot = CF*(g1qtot + g1qintregres + (1.5*cm.log(Q2/MuF2)-4.5-2*Zeta2)*(PDFpol(q,x,Q2)-PDFpol(-q,x,Q2)))

    def g1gint(z):
        g1gint = 1/z*0.5*((2*z-1)*(cm.log(Q2/MuF2*(1-z)/z)-1)+2*(1-z))*PDFpol(0,x/z,MuF2)
        return np.real(g1gint)

    g1gintres, g1ginterr = quad(g1gint, x, 1)

    return als(MuF2, 2)/(2*Pi)*(g1qtot + g1gintres)

def g5QCD(Q2, x, MuF2, q):

    def g5qintPP(z):
        g5qintPP = ((1/(1-z))*(2*cm.log(Q2/MuF2)-1.5) + 2*cm.log(1-z)/(1-z))
        g5qintPP = g5qintPP*((PDFpol(q,x/z,MuF2)-PDFpol(-q,x/z,MuF2))/z - PDFpol(q,x,MuF2) + PDFpol(-q,x,MuF2))
        return np.real(g5qintPP)

    g5qintPPres, g5qintPPerr = quad(g5qintPP, x, 1)

    def g5qintreg(z):
        g5qintreg = (-(1+z)*cm.log(Q2/MuF2*(1-z))-(1+z**2)/(1-z)*cm.log(z)+3+2*z)*(PDFpol(q,x/z,MuF2)-PDFpol(-q,x/z,MuF2))
        return np.real(g5qintreg)

    g5qintregres, g5qintregerr = quad(g5qintreg, x, 1)

    g5qtot = g5qintPPres - (PDFpol(q,x,Q2)-PDFpol(-q,x,Q2))*(0.5*(3-4*cm.log(Q2/MuF2))*cm.log(1-x)-2*li2(x))
    g5qtot = CF*(g5qtot + g5qintregres + (1.5*cm.log(Q2/MuF2)-4.5-2*Zeta2)*(PDFpol(q,x,Q2)-PDFpol(-q,x,Q2)))

    return als(MuF2, 2)/(2*Pi)*(g5qtot)

def gLQCD(Q2, x, MuF2, q):
    def gLqint(z):
        gLqint = 1/z*2*z*CF*(PDFpol(q,x/z,MuF2)-PDFpol(-q,x/z,MuF2))
        return np.real(gLqint)

    gLqintres, gLqinterr = quad(gLqint, x, 1)

    return gLqintres*als(MuF2, 2)/(2*Pi)

