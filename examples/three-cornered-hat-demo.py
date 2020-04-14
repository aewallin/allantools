"""
 Allan deviation tools
 Anders Wallin (anders.e.e.wallin "at" gmail.com)
 v1.0 2014 January

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Three-cornered-hat test
 
  See http://www.wriley.com/3-CornHat.htm
"""

import numpy
import matplotlib.pyplot as plt # only for plotting, not required for calculations

import allantools
from allantools import noise

def plotallan_phase(plt,y,rate,taus, style, label):
    (t2, ad, ade,adn) = allantools.oadev(y,data_type='phase',rate=rate,taus=taus)
    plt.loglog(t2, ad,  style, fillstyle='none', label=label,)

# plot a line with the slope alpha
def plotline(plt, alpha, amplitude, taus, style, labeltext):
    y = [ amplitude*pow(tt,alpha) for tt in taus]
    plt.loglog(taus,y,style, label=labeltext)
    
if __name__ == "__main__":
    print("allatools three-cornered-hat demo")
    # we test ADEV etc. by calculations on synthetic data
    # with known slopes of ADEV

    t = numpy.logspace( 0 ,4,50) # tau values from 1 to 1000
    plt.subplot(111, xscale="log", yscale="log")

    N=10000
    rate = 1.0
    # white phase noise => 1/tau ADEV
    
    ampl_A = 1.0 
    phaseA = ampl_A*numpy.random.randn(N)
    ampl_B = 5.0
    phaseB = ampl_B*numpy.random.randn(N)
    ampl_C = 10.0
    phaseC = ampl_C*numpy.random.randn(N)

    # measurements clk_I - clk_J
    phaseAB = phaseA - phaseB # [a-b for (a,b) in zip(phaseA,phaseB)]
    phaseBC = phaseB - phaseC # [b-c for (b,c) in zip(phaseB,phaseC)]
    phaseCA = phaseC - phaseA # [c-a for (c,a) in zip(phaseC,phaseA)]


    # theoretical ADEVs
    plotline(plt, -1.0, numpy.sqrt(3)*ampl_A, t, 'r--','A model')
    plotline(plt, -1.0, numpy.sqrt(3)*ampl_B, t, 'g--','B model')
    plotline(plt, -1.0, numpy.sqrt(3)*ampl_C, t, 'b--','C model')

    # simulated clocks
    plotallan_phase(plt, phaseA, 1, t, 'ro', 'true A phase')
    print("phaseA")
    plotallan_phase(plt, phaseB, 1, t, 'go', 'true B phase')
    print("phaseB")
    plotallan_phase(plt, phaseC, 1, t, 'bo', 'true C phase')
    print("phaseC")
    
    # clock differences
    plotallan_phase(plt, phaseAB, 1, t, 'r.', 'AB measurement')    
    print("phaseAB")
    plotallan_phase(plt, phaseBC, 1, t, 'g.', 'BC measurement')
    print("phaseBC")
    plotallan_phase(plt, phaseCA, 1, t, 'b.', 'CA measurement')
    print("phaseCA")
    
    # 3CH estimates
    (taus,devA,err_a,ns_ab) = allantools.three_cornered_hat_phase(phaseAB,phaseBC,phaseCA,rate,t, allantools.oadev)
    plt.loglog(taus, devA, 'rv', label='3-C-H estimate for A')
    
    (taus,devB,err_b,ns_ab) = allantools.three_cornered_hat_phase(phaseBC,phaseCA,phaseAB,rate,t, allantools.oadev)
    plt.loglog(taus, devB, 'gv', label='3-C-H estimate for B')
    
    (taus,devC,err_C,ns_ac) = allantools.three_cornered_hat_phase(phaseCA,phaseAB,phaseBC,rate,t, allantools.oadev)
    plt.loglog(taus, devC, 'bv', label='3-C-H estimate for C')
    
    print("TCH done.")
    plt.title('AllanTools three-cornered-hat example')
    plt.xlabel('Tau / s')
    plt.ylabel('OADEV')
    plt.legend()
    plt.grid()
    plt.show()
