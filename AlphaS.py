from scipy.interpolate import Rbf

######### Calculate the running of Alpha_S

def AlphaS(Mu2):
    if Mu2 < 4.18:
        Nfl = 4
    if 4.18 < Mu2 < 172.76:
        Nfl = 4
    if Mu2 > 172.76:
        Nfl = 6

    beta0 = 1/4*(11-2/3*Nfl)
    beta1 = 1/16*(102 - 38/3*Nfl)
    beta2 = 1/64*(2857/2-5033/18*Nfl+325/54*Nfl**2)

    return 1

####### 1/2/3/4-loop results for the running of alpha_S. Values generated with RunDec.m (arXiv:0004189)

####### Read out the numerical values of alpha_S for different values of the renormalization scale

alphasfile  = open('AlphaS1.txt', 'r')
readout = alphasfile.read().split()

alphasfile.close()

fulllist = [float(a) for a in readout]

Mupoints = fulllist[::2]
alphapoints = fulllist[1::2]

######## Interpolate the points to have a function alpha_S(Mu) for Mu up to 500GeV

def als(Mu, order):
    if order == 1:
        alsfit = Rbf(Mupoints[0:491], alphapoints[0:491])
        return alsfit(Mu)
    if order == 2:
        alsfit = Rbf(Mupoints[491:982], alphapoints[491:982])
        return alsfit(Mu)
    if order == 3:
        alsfit = Rbf(Mupoints[982:1473], alphapoints[982:1473])
        return alsfit(Mu)
    if order == 4:
        alsfit = Rbf(Mupoints[1473:1964], alphapoints[1473:1964])
        return alsfit(Mu)