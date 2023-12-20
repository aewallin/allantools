




    
import numpy
import matplotlib.pyplot as plt
import allantools as at

def Wpi(t, tau):
    # Pi counter
    if t > 0 and t <= tau:
        return 1.0/(tau)
    else:
        return 0.0

def Hpi(f, tau):
    # Pi counter frequency response
    return abs( numpy.sin( numpy.pi*f*tau ) / ( numpy.pi*f*tau) )

def Wallan(t, tau):
    # variance corresponding to Pi counter
    if t > 0 and t <= tau:
        return -1.0/(numpy.sqrt(2)*tau)
    elif t > tau and t <= 2*tau:
        return +1.0/(numpy.sqrt(2)*tau)
    else:
        return 0.0

def Wlambda(t, tau):
    # Lambda counter
    if t > 0 and t <= tau:
        return 1.0/(tau)*t
    elif t > tau and t <= 2*tau:
        return 1.0/(tau)-(1.0/tau)*(t-tau)
    else:
        return 0.0

def Hlambda(f, tau):
    return abs( pow( Hpi(f,tau), 2) )

    
def Wmod(t, tau):
    # variance corresponding to lambda counter
    # Dawkins (11)
    if t > 0 and t <= tau:
        return (-1.0/(numpy.sqrt(2)*pow(tau,2)))*t
    elif t > tau and t <= 2*tau:
        return (1.0/(numpy.sqrt(2)*pow(tau,2)))*(2*t-3*tau)
    elif t > 2*tau and t <= 3*tau:
        return (1.0/(numpy.sqrt(2)*pow(tau,2)))*(3*tau-t)
    else:
        return 0.0

def Womega(t, tau):
    # omega counter
    if t > 0 and t <= tau:
        return 3.0/(2*tau) - 6*pow(t-tau/2,2)/pow(tau,3)
    else:
        return 0.0

def WXomega(t, tau):
    # omega counter, phase weight
    if t > 0 and t <= tau:
        return -6.0/pow(tau,2) + 12*t/pow(tau,3)
    else:
        return 0.0


def Homega(f, tau):
    return abs( 3*numpy.sin( numpy.pi*f*tau ) /  pow( numpy.pi*f*tau,3) - 3*numpy.cos( numpy.pi*f*tau ) /  pow( numpy.pi*f*tau,2)   )
   
    
#%%
        
t = numpy.linspace(-5,8,50000)
f = numpy.linspace(0,5,500)

dt = min( numpy.diff(t) )
tau = 1.0
W_pi = [Wpi(x,tau) for x in t]
W_lam = [Wlambda(x,tau) for x in t]
W_om = [Womega(x,tau) for x in t]
WX_om = [WXomega(x,tau) for x in t]
#W_tri = [WT(x,tau) for x in t]
plt.figure()
plt.subplot(3,3,1)
plt.plot(t, W_pi,'b',label='$w_{\Pi}$')
plt.title('Weight for frequency data')
plt.grid()
plt.legend()
plt.xlim((-0.1,1.1))

plt.subplot(3,3,4)
plt.plot(t, W_lam,'m',label='$w_{\Lambda}$')
plt.xlim((-0.1,2.1))
plt.grid()
plt.legend()

plt.subplot(3,3,7)
plt.plot(t, W_om,'r',label='$w_{\Omega}$')
plt.xlim((-0.1,1.1))
plt.grid()
plt.legend()
plt.xlabel(' Time / t/$\\tau$')
plt.ylabel('Weight / 1/$\\tau$')

## phase weight
plt.subplot(3,3,2)
plt.arrow(0,0,0,1,width=0.02,length_includes_head=True)
plt.arrow(1,0,0,-1,width=0.02,length_includes_head=True)
plt.title('Weight for phase data')
plt.grid()
plt.legend()
plt.xlim((-0.1,1.1))
plt.ylim((-1.1,1.1))

plt.subplot(3,3,5)
plt.plot(t[1:], numpy.diff(W_lam)/min(numpy.diff(t)),'m',label='$w_{\Lambda}^x$')
plt.xlim((-0.1,2.1))
plt.grid()
plt.legend()

plt.subplot(3,3,8)
plt.plot(t, WX_om,'r',label='$w_{\Omega}^x$')
plt.xlim((-0.1,1.1))
plt.grid()
plt.legend()
plt.xlabel(' Time / t/$\\tau$')
plt.ylabel('Weight / 1/$\\tau$')



# frequency response
plt.subplot(3,3,3)
plt.plot(f, Hpi(f,1.0),'b',label='$abs( FFT( w_{\Pi} ) )$')
plt.title('Frequency response of counter')
plt.grid()
plt.legend()
plt.ylim((0,1.05))
#plt.xlim((-0.1,1.1))
    
plt.subplot(3,3,6)
plt.plot(f, Hlambda(f,1.0),'m',label='$abs( FFT( w_{\Lambda} ) )$')
plt.ylim((0,1.05))
plt.grid()
plt.legend()

plt.subplot(3,3,9)
plt.plot(f, Homega(f,1.0),'r',label='$abs( FFT( w_{\Omega} ) )$')
#plt.xlim((-0.1,1.1))
plt.grid()
plt.legend()
plt.ylim((0,1.05))
#plt.xlabel(' Time / t/$\\tau$')
#plt.ylabel('Weight / 1/$\\tau$')
plt.xlabel(' Normalized frequency / $f$ $\\tau$')

plt.show()
