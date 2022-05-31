from Counterterms import *

#### Z vertex functions - normalized to the vector/axial born coupling

def gVeNLO(p12,MuR2):
    gVe1 = (4*(-1 + 2*CW2)*(MW2 + 2*p12)*A0fin(MW2,MuR2) -
            2*(MZ2 + 2*p12)*(-1 + 6*SW2 - 12*SW2**2 + 16*SW2**4)*A0fin(MZ2,MuR2) +
            MZ2*(4*CW2*p12 + 2*p12*(-1 + 6*SW2 - 12*SW2**2 + 16*SW2**4) +
            (2*CW2*(2*MW2 + 3*p12) + (2*MZ2 + 3*p12)*(-1 + 6*SW2 - 12*SW2**2 + 16*SW2**4))*
            B0fin(p12, 0, 0,MuR2)) - 4*CW2**2*MZ2*(2*MW2 + p12)*B0fin(p12, MW2, MW2, MuR2) +
            4*CW2*MZ2*(MW2**2 + 2*MW2*p12 + p12**2)*C01(p12, MW2) +
            2*MZ2*(MZ2**2 + 2*MZ2*p12 + p12**2)*(-1 + 6*SW2 - 12*SW2**2 + 16*SW2**4)*C01(p12, MZ2) +
            8*CW2**2*MW2*MZ2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (gVe1*aem/(32*MW2*p12*Pi*SW2) + (dZZZ1(MuR2) - (2*dSW1(MuR2))/SW - (6*dSW1(MuR2)*SW)/CW2 - 4*CW*dZAZ1(MuR2)*SW
            - 4*dZZZ1(MuR2)*SW2 +dZe1(MuR2)*(2 - 8*SW2) + (2- 4*SW2)*dZeL(MuR2) - 4*SW2*dZeR(MuR2))/4)

def gAeNLO(p12,MuR2):
    gAe1 = (4*(-1 + 2*CW2)*(MW2 + 2*p12)*A0fin(MW2,MuR2) + 2*(MZ2 + 2*p12)*(1 - 6*SW2 + 12*SW2**2)*
            A0fin(MZ2,MuR2) + MZ2*(4*CW2*p12 - 2*p12*(1 - 6*SW2 + 12*SW2**2) +
            (2*CW2*(2*MW2 + 3*p12) + (-2*MZ2 - 3*p12)*(1 - 6*SW2 + 12*SW2**2))*B0fin(p12, 0, 0,MuR2)) -
            4*CW2**2*MZ2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) + 4*(MW2**3 + 2*MW2**2*p12 + MW2*p12**2)*
            C01(p12, MW2) - 2*MZ2*(MZ2**2 + 2*MZ2*p12 + p12**2)*(1 - 6*SW2 + 12*SW2**2)*
            C01(p12, MZ2) + 8*CW2*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (gAe1*aem/(32*MW2*p12*Pi*SW2) + (2*dZe1(MuR2) + dZZZ1(MuR2) + 2*(-(dSW1(MuR2)/SW) + (dSW1(MuR2)*SW)/CW2
            + (1-2*SW2)*dZeL(MuR2) + 2*SW2*dZeR(MuR2)))/4)

def gVuNLO(p12,MuR2):
    gVu1 = (2*MZ2*p12*(27 - 108*SW2 + 144*SW2**2 - 128*SW2**3 + 18*CW2*(-3 + 2*SW2)) -
            36*(MW2 + 2*p12)*(-3 + 6*CW2 + 2*SW2)*A0fin(MW2,MuR2) +
            2*(MZ2 + 2*p12)*(-27 + 108*SW2 - 144*SW2**2 + 128*SW2**3)*A0fin(MZ2,MuR2) +
            MZ2*(18*CW2*(2*MW2 + 3*p12)*(-3 + 2*SW2) - (2*MZ2 + 3*p12)*
            (-27 + 108*SW2 - 144*SW2**2 + 128*SW2**3))*B0fin(p12, 0, 0,MuR2) +
            108*CW2*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) + 36*MW2*(MW2 + p12)**2*(-3 + 2*SW2)*
            C01(p12, MW2) - 2*MZ2*(MZ2 + p12)**2*(-27 + 108*SW2 - 144*SW2**2 +
            128*SW2**3)*C01(p12, MZ2) - 216*CW2*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (gVu1*aem/(864*MW2*p12*Pi*SW2) + ((6*dSW1(MuR2))/SW + (10*dSW1(MuR2)*SW)/CW2 + 8*CW*dZAZ1(MuR2)*SW
            + 2*dZe1(MuR2)*(-3 + 8*SW2) +dZZZ1(MuR2)*(-3 + 8*SW2) + (8*SW2- 6)*dZuL(MuR2) + 8*SW2*dZuR(MuR2))/12)

def gAuNLO(p12,MuR2):
    gAu1 = (2*MZ2*p12*(3 - 12*SW2 + 16*SW2**2 + CW2*(-6 + 4*SW2)) -
            4*(MW2 + 2*p12)*(-3 + 6*CW2 + 2*SW2)*A0fin(MW2,MuR2) -
            2*(MZ2 + 2*p12)*(3 - 12*SW2 + 16*SW2**2)*A0fin(MZ2,MuR2) +
            MZ2*(2*CW2*(2*MW2 + 3*p12)*(-3 + 2*SW2) + (2*MZ2 + 3*p12)*(3 - 12*SW2 + 16*SW2**2))*
            B0fin(p12, 0, 0,MuR2) + 12*CW2*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) +
            4*MW2*(MW2 + p12)**2*(-3 + 2*SW2)*C01(p12, MW2) +
            2*MZ2*(MZ2 + p12)**2*(3 - 12*SW2 + 16*SW2**2)*C01(p12, MZ2) -
            24*CW2*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (gAu1*aem/(96*MW2*p12*Pi*SW2) + (-6*dZe1(MuR2) - 3*dZZZ1(MuR2) + (6*dSW1(MuR2))/SW
            - (6*dSW1(MuR2)*SW)/CW2 + (8*SW2- 6)*dZuL(MuR2) - 8*SW2*dZuR(MuR2))/12)

def gVdNLO(p12,MuR2):
    gVd1 = (36*(MW2 + 2*p12)*(-3 + 6*CW2 + 4*SW2)*A0fin(MW2,MuR2) -
            2*(MZ2 + 2*p12)*(-27 + 54*SW2 - 36*SW2**2 + 16*SW2**3)*A0fin(MZ2,MuR2) +
            MZ2*(-3 + 4*SW2)*(2*p12*(9 - 18*CW2 - 6*SW2 + 4*SW2**2) +
            (-18*CW2*(2*MW2 + 3*p12) + (2*MZ2 + 3*p12)*(9 - 6*SW2 + 4*SW2**2))*B0fin(p12, 0, 0,MuR2)) -
            108*CW2*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) - 36*MW2*(MW2 + p12)**2*(-3 + 4*SW2)*
            C01(p12, MW2) + 2*MZ2*(MZ2 + p12)**2*(-27 + 54*SW2 - 36*SW2**2 + 16*SW2**3)*
            C01(p12, MZ2) + 216*CW2*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (gVd1*aem/(864*MW2*p12*Pi*SW2) + (-(6*CW2*dSW1(MuR2) + 2*dSW1(MuR2)*SW2 + 4*CW**3*dZAZ1(MuR2)*SW2 +
            CW2*SW*(dZZZ1(MuR2)*(-3 + 4*SW2) + dZe1(MuR2)*(-6 + 8*SW2) + (4*SW2- 6)*dZdL(MuR2) + 4*SW2*dZdR(MuR2)))/(12*CW2*SW)))

def gAdNLO(p12,MuR2):
    gAd1 = (-2*MZ2*p12*(3 - 6*SW2 + 4*SW2**2 + CW2*(-6 + 8*SW2)) +
            4*(MW2 + 2*p12)*(-3 + 6*CW2 + 4*SW2)*A0fin(MW2,MuR2) + 2*(MZ2 + 2*p12)*(3 - 6*SW2 + 4*SW2**2)*
            A0fin(MZ2,MuR2) - MZ2*(2*CW2*(2*MW2 + 3*p12)*(-3 + 4*SW2) +
            (2*MZ2 + 3*p12)*(3 - 6*SW2 + 4*SW2**2))*B0fin(p12, 0, 0,MuR2) -
            12*CW2*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) - 4*MW2*(MW2 + p12)**2*(-3 + 4*SW2)*
            C01(p12, MW2) - 2*MZ2*(MZ2 + p12)**2*(3 - 6*SW2 + 4*SW2**2)*C01(p12, MZ2) +
            24*CW2*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (gAd1*aem/(96*MW2*p12*Pi*SW2) + (-6*CW2*dSW1(MuR2) + 6*dSW1(MuR2)*SW2 + CW2*SW*(6*dZe1(MuR2) + 3*dZZZ1(MuR2)
            + (6-4*SW2)*dZdL(MuR2) + 4*SW2*dZdR(MuR2)))/(12*CW2*SW))

#### Photon vertex functions - normalized to the vector born coupling

def hVeNLO(p12,MuR2):
    hVe1 = (-4*(MW2 + 2*p12)*A0fin(MW2,MuR2) - 2*(MZ2 + 2*p12)*(1 - 4*SW2 + 8*SW2**2)*A0fin(MZ2,MuR2) +
            MZ2*(1 - 4*SW2 + 8*SW2**2)*(2*p12 + (2*MZ2 + 3*p12)*B0fin(p12, 0, 0,MuR2)) +
            2*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) + 2*MZ2*(MZ2 + p12)**2*(1 - 4*SW2 + 8*SW2**2)*
            C01(p12, MZ2) - 4*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (hVe1*aem/(32*MW2*p12*Pi*SW2) + (dZZA1(MuR2) - 4*dZZA1(MuR2)*SW2 - 4*CW*SW*(dZAA1(MuR2) + 2*dZe1(MuR2)
            + dZeL(MuR2) +dZeR(MuR2)))/(8*CW*SW))

def hAeNLO(p12,MuR2):
    hAe1 = (-4*(MW2 + 2*p12)*A0fin(MW2,MuR2) + 2*(MZ2 + 2*p12)*(-1 + 4*SW2)*A0fin(MZ2,MuR2) +
            MZ2*(-1 + 4*SW2)*(-2*p12 - (2*MZ2 + 3*p12)*B0fin(p12, 0, 0,MuR2)) +
            2*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) - 2*MZ2*(MZ2 + p12)**2*(-1 + 4*SW2)*
            C01(p12, MZ2) - 4*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (hAe1*aem/(32*MW2*p12*Pi*SW2) + (dZZA1(MuR2) + 4*CW*SW*(-dZeL(MuR2) + dZeR(MuR2)))/(8*CW*SW))

def hVuNLO(p12,MuR2):
    hVu1 = (2*MZ2*p12*(-9 + 9*CW2 + 24*SW2 - 32*SW2**2) + 36*(MW2 + 2*p12)*A0fin(MW2,MuR2) +
            2*(MZ2 + 2*p12)*(9 - 24*SW2 + 32*SW2**2)*A0fin(MZ2,MuR2) +
            MZ2*(9*CW2*(2*MW2 + 3*p12) - (2*MZ2 + 3*p12)*(9 - 24*SW2 + 32*SW2**2))*
            B0fin(p12, 0, 0,MuR2) - 27*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) +
            18*MW2*(MW2 + p12)**2*C01(p12, MW2) - 2*MZ2*(MZ2 + p12)**2*(9 - 24*SW2 + 32*SW2**2)*
            C01(p12, MZ2) + 54*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (hVu1*aem/(432*MW2*p12*Pi*SW2) + ((dZZA1(MuR2)*(-3 + 8*SW2))/(CW*SW) + 8*(dZAA1(MuR2) + 2*dZe1(MuR2)
            + dZuL(MuR2) + dZuR(MuR2)))/24)

def hAuNLO(p12,MuR2):
    hAu1 = (2*MZ2*p12*(-3 + 3*CW2 + 8*SW2) + 12*(MW2 + 2*p12)*A0fin(MW2,MuR2) -
            2*(MZ2 + 2*p12)*(-3 + 8*SW2)*A0fin(MZ2,MuR2) +
            MZ2*(CW2*(6*MW2 + 9*p12) + (2*MZ2 + 3*p12)*(-3 + 8*SW2))*B0fin(p12, 0, 0,MuR2) -
            9*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) + 6*MW2*(MW2 + p12)**2*C01(p12, MW2) +
            2*MZ2*(MZ2 + p12)**2*(-3 + 8*SW2)*C01(p12, MZ2) +
            18*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (hAu1*aem/(144*MW2*p12*Pi*SW2) + ((-3*dZZA1(MuR2))/(CW*SW) + 8*dZuL(MuR2) - 8*dZuR(MuR2))/24)

def hVdNLO(p12,MuR2):
    hVd1 = (2*MZ2*p12*(9 - 36*CW2 - 12*SW2 + 8*SW2**2) - 36*(MW2 + 2*p12)*A0fin(MW2,MuR2) -
            2*(MZ2 + 2*p12)*(9 - 12*SW2 + 8*SW2**2)*A0fin(MZ2,MuR2) +
            MZ2*(-36*CW2*(2*MW2 + 3*p12) + (2*MZ2 + 3*p12)*(9 - 12*SW2 + 8*SW2**2))*
            B0fin(p12, 0, 0,MuR2) + 54*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) -
            72*MW2*(MW2 + p12)**2*C01(p12, MW2) + 2*MZ2*(MZ2 + p12)**2*(9 - 12*SW2 + 8*SW2**2)*
            C01(p12, MZ2) - 108*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (hVd1*aem/(864*MW2*p12*Pi*SW2) + (dZZA1(MuR2)*(3 - 4*SW2) - 4*CW*SW*(dZAA1(MuR2) + 2*dZe1(MuR2)
            + dZdL(MuR2) + dZdR(MuR2)))/(24*CW*SW))

def hAdNLO(p12,MuR2):
    hAd1 = (-2*MZ2*p12*(-3 + 12*CW2 + 4*SW2) - 12*(MW2 + 2*p12)*A0fin(MW2,MuR2) +
            2*(MZ2 + 2*p12)*(-3 + 4*SW2)*A0fin(MZ2,MuR2) -
            MZ2*(12*CW2*(2*MW2 + 3*p12) + (2*MZ2 + 3*p12)*(-3 + 4*SW2))*B0fin(p12, 0, 0,MuR2) +
            18*MW2*(2*MW2 + p12)*B0fin(p12, MW2, MW2,MuR2) - 24*MW2*(MW2 + p12)**2*C01(p12, MW2) -
            2*MZ2*(MZ2 + p12)**2*(-3 + 4*SW2)*C01(p12, MZ2) -
            36*MW2**2*(MW2 + 2*p12)*C02(p12, MW2, MW2))
    return (hAd1*aem/(288*MW2*p12*Pi*SW2) + ((3*dZZA1(MuR2))/(CW*SW) - 4*dZdL(MuR2) + 4*dZdR(MuR2))/24)


#### Actual Vertex functions for the Z vertex

def gVe(Q2, MuR2, n):
    if n == 0:
        return -0.5+2*SW2
    else:
        return gVeNLO(Q2,MuR2)

def gAe(Q2, MuR2, n):
    if n == 0:
        return -0.5
    else:
        return gAeNLO(Q2,MuR2)

def gVq(Q2, MuR2, q, n):
    if n == 0:
        if (q % 2) == 0:
            return 0.5 - 4*SW2/3
        else:
            return -0.5 + 2*SW2/3
    else:
        if (q % 2) == 0:
            return gVuNLO(Q2,MuR2)
        else:
            return gVdNLO(Q2,MuR2)

def gAq(Q2, MuR2, q, n):
    if n == 0:
        if (q % 2) == 0:
            return 0.5
        else:
            return -0.5
    else:
        if (q % 2) == 0:
            return gAuNLO(Q2,MuR2)
        else:
            return gAdNLO(Q2,MuR2)

#### Actual Vertex functions for the Photon vertex

def hVe(Q2, MuR2, n):
    if n == 0:
        return -1
    else:
        return hVeNLO(Q2,MuR2)

def hAe(Q2, MuR2, n):
    if n == 0:
        return 0
    else:
        return hAeNLO(Q2,MuR2)

def hVq(Q2, MuR2, q, n):
    if n == 0:
        if (q % 2) == 0:
            return 2/3
        else:
            return -1/3
    else:
        if (q % 2) == 0:
            return hVuNLO(Q2,MuR2)
        else:
            return hVdNLO(Q2,MuR2)

def hAq(Q2, MuR2, q, n):
    if n == 0:
        return 0
    else:
        if (q % 2) == 0:
            return hAuNLO(Q2,MuR2)
        else:
            return hAdNLO(Q2,MuR2)

# W Vertex functions

def kenu(p12,MuR2):
    kenu = (4*(MW2 + 2*p12)*A0fin(MW2, MuR2) + 2*(MZ2 + 2*p12)*(1 - 2*SW2 + 2*SW2**2)*
                  A0fin(MZ2, MuR2) - MZ2*(2*MZ2 + 3*p12)*(-1 + 2*SW2)*B0fin(p12, 0, 0, MuR2) -
                  4*MW2*(MW2 + p12)*SW2*B0fin(p12, 0, MW2, MuR2) - 4*CW2*MW2*(MW2 + MZ2 + p12)*
                  B0fin(p12, MW2, MZ2, MuR2) + 8*MW2**2*p12*SW2*C01(p12, MW2) -
                  2*MZ2*(MZ2 + p12)**2*(-1 + 2*SW2)*C01(p12, MZ2) +
                  8*CW2*MW2*(MZ2*p12 + MW2*(MZ2 + p12))*C02(p12, MW2, MZ2))

    return kenu*aem/(8*MW2*Pi*p12*SW2)

def kdu(p12,MuR2):
    kdu = (36*(MW2 + 2*p12)*A0fin(MW2, MuR2) + 2*(MZ2 + 2*p12)*(9 - 18*SW2 + 10*SW2**2)*
                 A0fin(MZ2, MuR2) + 40*MW2*p12*SW2*B0fin(0, 0, 0, MuR2) +
                 MZ2*(3*p12*(9 - 10*SW2) + 2*MZ2*(9 - 18*SW2 + 8*SW2**2))*B0fin(p12, 0, 0, MuR2) -
                 36*MW2*(MW2 + p12)*SW2*B0fin(p12, 0, MW2, MuR2) - 36*CW2*MW2*(MW2 + MZ2 + p12)*
                 B0fin(p12, MW2, MZ2, MuR2) + 16*MW2*p12**2*SW2*C0PV[0, 0, p12, 0, 0, 0] +
                 72*MW2**2*p12*SW2*C01(p12, MW2) + 2*MZ2*(MZ2 + p12)**2*
                 (9 - 18*SW2 + 8*SW2**2)*C01(p12, MZ2) +
                 72*CW2*MW2*(MZ2*p12 + MW2*(MZ2 + p12))*C02(p12, MW2, MZ2))

    return kdu*aem/(72*MW2*p12*Pi*SW2)