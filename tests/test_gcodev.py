"""
  Gcodev test for allantools (https://github.com/aewallin/allantools)

"""
import math
# import time
# import sys
import pytest
import numpy
import allantools


# SigmaTheta v3.0-29-g053a Wed Oct 11 12:59:22 2023 +0000
#
# ./GCoDev data/freqCA.txt data/freqAB.txt 
# cannot open config file /usr/share/sigmatheta/SigmaTheta.conf
# using hardcoded defaults.
# AF                          GCoDev( freqCA, freqAB )
SigmaTheta_A = [ ( 1.0000000000000000e+00, -1.5822698768801253e+00),
                 ( 2.0000000000000000e+00, -7.1516620564809252e-01),
                 ( 4.0000000000000000e+00, -3.6211155665136197e-01),
                 ( 8.0000000000000000e+00, -1.7060913297684019e-01),
                 ( 1.6000000000000000e+01, -1.0714277882857023e-01),
                 ( 3.2000000000000000e+01, -5.0828253381873023e-02),
                 ( 6.4000000000000000e+01, -2.3379022940577421e-02),
                 ( 1.2800000000000000e+02, -1.3214052398385969e-02),
                 ( 2.5600000000000000e+02, -5.9325361833646468e-03),
                 ( 5.1200000000000000e+02, -3.1859793758841873e-03),
                 ( 1.0240000000000000e+03, -1.1872594305012082e-03),
                 ( 2.0480000000000000e+03, -7.9181618114116986e-04),
                 ( 4.0960000000000000e+03, -3.1864022889666255e-04),
                 ( 8.1920000000000000e+03, -1.8946156459394094e-04),
                 ( 1.6384000000000000e+04, -9.2401802887092055e-05),
                 ( 3.2768000000000000e+04, -4.7453355340743471e-05)]

# SigmaTheta v3.0-29-g053a Wed Oct 11 12:59:22 2023 +0000
#
# ./GCoDev data/freqBC.txt data/freqAB.txt 
# cannot open config file /usr/share/sigmatheta/SigmaTheta.conf
# using hardcoded defaults.
# AF                          GCoDev( freqBC, freqCA )
SigmaTheta_B = [ (   1.0000000000000000e+00, -8.6572265057967588e+00),
                (  2.0000000000000000e+00, -4.3432005973471712e+00),
                (  4.0000000000000000e+00, -2.1721589652062332e+00),
                (  8.0000000000000000e+00, -1.0880738273515811e+00),
                (  1.6000000000000000e+01, -5.4146914222235154e-01),
                (  3.2000000000000000e+01, -2.7216329996649430e-01),
                (  6.4000000000000000e+01, -1.3536179632154205e-01),
                (  1.2800000000000000e+02, -6.7521150727460719e-02),
                (  2.5600000000000000e+02, -3.3923783921641061e-02),
                (  5.1200000000000000e+02, -1.6903085282002191e-02),
                (  1.0240000000000000e+03, -8.5454849855555823e-03),
                (  2.0480000000000000e+03, -4.2376272243631599e-03),
                (  4.0960000000000000e+03, -2.1168485094171896e-03),
                (  8.1920000000000000e+03, -1.0600373993369303e-03),
                (  1.6384000000000000e+04, -5.2776620159470755e-04),
                (  3.2768000000000000e+04, -2.6352955508714208e-04)]

# SigmaTheta v3.0-29-g053a Wed Oct 11 12:59:22 2023 +0000
# 
# /GCoDev data/freqBC.txt data/freqCA.txt 
# cannot open config file /usr/share/sigmatheta/SigmaTheta.conf
# using hardcoded defaults.
# AF                          GCoDev( freqBC, freqCA )
SigmaTheta_C = [ ( 1.0000000000000000e+00, -1.7326873626967341e+01),
                (  2.0000000000000000e+00, -8.6728637511045310e+00),
                (  4.0000000000000000e+00, -4.3259160453426837e+00),
                (  8.0000000000000000e+00, -2.1672253715364658e+00),
                (  1.6000000000000000e+01, -1.0832324591986147e+00),
                (  3.2000000000000000e+01, -5.4196999788713551e-01),
                (  6.4000000000000000e+01, -2.7112486792251950e-01),
                (  1.2800000000000000e+02, -1.3560263809161266e-01),
                (  2.5600000000000000e+02, -6.7501039944739105e-02),
                (  5.1200000000000000e+02, -3.3778774118065999e-02),
                (  1.0240000000000000e+03, -1.6967984008809077e-02),
                (  2.0480000000000000e+03, -8.4465947630849424e-03),
                (  4.0960000000000000e+03, -4.2428213734232445e-03),
                (  8.1920000000000000e+03, -2.1174831713595926e-03),
                (  1.6384000000000000e+04, -1.0558534636972263e-03),
                (  3.2768000000000000e+04, -5.2745664531139290e-04)]
                  

class TestGcodev():
    def generate_data(self):
        N = 100000
        rate = 1.0
        # white phase noise => 1/tau ADEV
        # these are the 'true' phases of the oscillators, which are not observable.
        # we can only measure phases between two oscillators.
        numpy.random.seed(42)
        ampl_A = 1.0
        phaseA = ampl_A*numpy.random.randn(N)
        ampl_B = 5.0
        phaseB = ampl_B*numpy.random.randn(N)
        ampl_C = 10.0
        phaseC = ampl_C*numpy.random.randn(N)
        print('Generated random data:')
        print(phaseA[-1],phaseB[-1],phaseC[-1])
        # seed = 42 output:
        # 0.12006294082414522 -5.324739839420351 -3.578461741191343
        
        # measurements clk_I - clk_J
        phaseAB = phaseA - phaseB  # [a-b for (a,b) in zip(phaseA,phaseB)]
        phaseBC = phaseB - phaseC  # [b-c for (b,c) in zip(phaseB,phaseC)]
        phaseCA = phaseC - phaseA  # [c-a for (c,a) in zip(phaseC,phaseA)]
        return phaseAB, phaseBC, phaseCA
    
    def check_devs(self, dev2, dev1, soft=False):
        rel_error = (dev2-dev1)/dev1
        tol = 1e-4
        verbose = 1
    
        if (abs(rel_error) < tol):
            if verbose:
                print(" OK   %0.6f \t    %0.6f \t %0.6f" % (dev1, dev2, rel_error))
            return True
        else:
            print(" ERROR   %0.6f \t %0.6f \t %0.6f" % (dev1, dev2, rel_error))
            print(" bias= %.4f" % pow(dev2/dev1, 2))
            if soft:
                return True
            return False
        
    def test_gcodev(self):
        phaseAB, phaseBC, phaseCA = self.generate_data()
        # GCODEV estimates
        (tauA, devA, err_a, ns_ab) = allantools.gcodev(phaseAB, phaseCA, rate=1, taus='octave')
        print('Gcodev A', devA)
        for (st, A) in zip (SigmaTheta_A, devA):
            self.check_devs(A, -st[1])
        
        (tauB, devB, err_b, ns_b) = allantools.gcodev(phaseAB, phaseBC, rate=1, taus='octave')
        print('Gcodev B', devB)
        for (st, B) in zip (SigmaTheta_B, devB):
            self.check_devs(B, -st[1])
            
        (tauC, devC, err_c, ns_c) = allantools.gcodev(phaseCA, phaseBC, rate=1, taus='octave')
        print('Gcodev C', devC)
        for (st, C) in zip (SigmaTheta_C, devC):
            self.check_devs(C, -st[1])
            
        print("Gcodev done.")
    
        





if __name__ == "__main__":
    t = TestGcodev()
    t.test_gcodev()
    #pytest.main(["test_gcodev.py"])
