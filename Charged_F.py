from Counterterms import *

#Charged current structure functions

# Lepton Vertex Corrections
def F1euLEP(p12, MuR2):
    F1eu = (((8*MW2*(CW2*(MZ2 + p12) + p12*(-2 + SW2)) -
                2*MZ2*(4*MZ2*(-1 + SW2)**2 + p12*(7 - 14*SW2 + 8*SW2**2)) +
                16*MW2**2*SW2*cm.log(1 - p12/MW2) -
                2*(MZ2*(-1 + 2*SW2)*(2*p12 + (2*MZ2 + 3*p12)*B0fin(p12, 0, 0,MuR2)) +
                4*MW2*(MW2 + p12)*SW2*B0fin(p12, 0, MW2,MuR2) + 4*CW2*MW2*(MW2 + MZ2 + p12)*
                B0fin(p12, MW2, MZ2,MuR2) + 2*MZ2*(MZ2 + p12)**2*(-1 + 2*SW2)*
                C01(p12, MZ2) - 8*CW2*MW2*(MZ2*p12 + MW2*(MZ2 + p12))*
                C02(p12, MW2, MZ2) + 4*MW2*(MW2 + 2*p12)*(-1 + cm.log(MW2/MuR2)) +
                2*MZ2*(MZ2 + 2*p12)*(1 - 2*SW2 + 2*SW2**2)*(-1 + cm.log(MZ2/MuR2)) -
                8*MW2**2*SW2*(cm.log(1 - p12/MW2)*(cm.log(MW2/MuR2) + cm.log(1 - p12/MW2)) + li2(p12/MW2))))))

    F1eu = F1eu*aem/(16*MW2*p12*Pi*SW2)

    return F1eu

def FLeuLEP(p12, MuR2):
    Fleu = (-4*CW2**2*MZ2**2 + 4*CW2*MW2*(MZ2 + p12) + 4*MW2*p12*(-2 + SW2) +
            MZ2*p12*(-7 + 14*SW2 - 8*SW2**2) + 8*MW2**2*SW2*cm.log(1 - p12/MW2))
    Fleu = -Fleu*aem/(4*MW2*p12*Pi*SW2)

    return Fleu

def F3euLEP(i, p12, MuR2):
    F3eu = (-(MZ2*(2*MZ2 + 3*p12)*(-1 + 2*SW2)*B0fin(p12, 0, 0,MuR2)) -
                 4*MW2*(MW2 + p12)*SW2*B0fin(p12, 0, MW2,MuR2) - 4*CW2*MW2*(MW2 + MZ2 + p12)*
                 B0fin(p12, MW2, MZ2,MuR2) - 2*MZ2*(MZ2 + p12)**2*(-1 + 2*SW2)*
                 C01(p12, MZ2) + 8*CW2*MW2*(MZ2*p12 + MW2*(MZ2 + p12))*
                 C02(p12, MW2, MZ2) - 4*MW2*(MW2 + 2*p12)*(-1 + cm.log(MW2/MuR2)) -
                 2*MZ2*(MZ2 + 2*p12)*(1 - 2*SW2 + 2*SW2**2)*(-1 + cm.log(MZ2/MuR2)) +
                 8*MW2**2*SW2*(cm.log(1 - p12/MW2)*(cm.log(MW2/MuR2) + cm.log(1 - p12/MW2)) + li2(p12/MW2)))
    F3eu = F3eu*aem/(4*MW2*p12*Pi*SW2)

    if (abs(i) % 2) == 0:
        return (F3eu*LepType + (1-LepType)/2*aem/(2*Pi)*(1+2*SW2)/SW2)*abs((LepType+i/abs(i))/2)
    elif (abs(i) != 2):
        return -F3eu*LepType*abs((LepType-i/abs(i))/2)



# Quark Vertex Corrections
def F1euQuark(p12, MuR2):

    F1eu = (72*CW2*MW2*(MW2 + MZ2 + p12) - 72*MW2*(MW2 + 2*p12) +
                72*MW2*(MW2 + p12)*SW2 + 6*MZ2*p12*(-9 + 10*SW2) -
                4*MZ2**2*(9 - 18*SW2 + 8*SW2**2) - 4*MZ2*(MZ2 + 2*p12)*(9 - 18*SW2 + 10*SW2**2) +
                32*MW2*p12*SW2*cm.log(p12/MuR2) + 144*MW2**2*SW2*cm.log(1 - p12/MW2) -
                2*(MZ2*(2*p12*(-9 + 10*SW2) + (3*p12*(-9 + 10*SW2) -
                2*MZ2*(9 - 18*SW2 + 8*SW2**2))*B0fin(p12, 0, 0,MuR2)) +
                36*MW2*(MW2 + p12)*SW2*B0fin(p12, 0, MW2,MuR2) + 36*CW2*MW2*(MW2 + MZ2 + p12)*
                B0fin(p12, MW2, MZ2,MuR2) - 2*MZ2*(MZ2 + p12)**2*(9 - 18*SW2 + 8*SW2**2)*
                C01(p12, MZ2) - 72*CW2*MW2*(MZ2*p12 + MW2*(MZ2 + p12))*
                C02(p12, MW2, MZ2) + 36*MW2*(MW2 + 2*p12)*(-1 + cm.log(MW2/MuR2)) +
                2*MZ2*(MZ2 + 2*p12)*(9 - 18*SW2 + 10*SW2**2)*(-1 + cm.log(MZ2/MuR2)) +
                (4*MW2*p12*SW2*(7*Pi**2 - 6*cm.log(p12/MuR2)**2))/3 -
                72*MW2**2*SW2*(cm.log(1 - p12/MW2)*(cm.log(MW2/MuR2) + cm.log(1 - p12/MW2)) + li2(p12/MW2))))

    F1eu = F1eu*aem/(144*MW2*p12*Pi*SW2)

    return F1eu


def FLeuQuark(p12, MuR2):
    Fleu = -((36*MW2*(CW2*(MZ2 + p12) + p12*(-2 + SW2)) -
                   MZ2*(36*CW2**2*MZ2 + p12*(63 - 102*SW2 + 40*SW2**2)) +
                   8*MW2*SW2*(2*p12*cm.log(p12/MuR2) + 9*MW2*cm.log(1 - p12/MW2))))
    Fleu = -Fleu*aem/(36*MW2*p12*Pi*SW2)

    return Fleu

def F3euQuark(i, p12, MuR2):
    F3eu = ((MZ2*(3*p12*(9 - 10*SW2) + 2*MZ2*(9 - 18*SW2 + 8*SW2**2))*B0fin(p12, 0, 0,MuR2) -
                 36*MW2*(MW2 + p12)*SW2*B0fin(p12, 0, MW2,MuR2) - 36*CW2*MW2*(MW2 + MZ2 + p12)*
                 B0fin(p12, MW2, MZ2,MuR2) + 2*MZ2*(MZ2 + p12)**2*(9 - 18*SW2 + 8*SW2**2)*
                 C01(p12, MZ2) + 72*CW2*MW2*(MZ2*p12 + MW2*(MZ2 + p12))*
                 C02(p12, MW2, MZ2) - 36*MW2*(MW2 + 2*p12)*(-1 + cm.log(MW2/MuR2)) -
                 2*MZ2*(MZ2 + 2*p12)*(9 - 18*SW2 + 10*SW2**2)*(-1 + cm.log(MZ2/MuR2)) +
                 (4*MW2*p12*SW2*(-7*Pi**2 + 6*cm.log(p12/MuR2)**2))/3 +
                 72*MW2**2*SW2*(cm.log(1 - p12/MW2)*(cm.log(MW2/MuR2) + cm.log(1 - p12/MW2)) + li2(p12/MW2))))
    F3eu = F3eu*aem/(36*MW2*p12*Pi*SW2)

    if (abs(i) % 2) == 0:
            return (F3eu*LepType + (1-LepType)/2*7*aem/(9*Pi))*abs((LepType+i/abs(i))/2)
    elif (abs(i) != 2):
        return (-F3eu*LepType+ (1+LepType)/2*7*aem/(9*Pi))*abs((LepType-i/abs(i))/2)


