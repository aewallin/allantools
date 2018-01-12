import allantools as at
from matplotlib import pyplot as plt
import numpy as np

def b1_noise_id(x, af, rate):
    """ R(n) noise identification algorithm
    
        based on expected ratio of MVAR/AVAR
    """
    #print "len(x) ", len(x), "oadev",af*rate
    (taus,devs,errs,ns) = at.adev(x,taus=[af*rate],data_type="phase", rate=rate) 
    #print devs
    oadev_x = devs[0]
    y = np.diff(x)
    #y_dec = y[0:len(y):af]
    y_cut = np.array( y[:len(y)-(len(y)%af)] ) # cut to length
    assert len(y_cut)%af == 0
    y_shaped = y_cut.reshape( ( int(len(y_cut)/af), af) )
    y_averaged = np.average(y_shaped,axis=1) # average
        
    #x_dec = x[0:len(x):af]
    var = np.var(y_averaged, ddof=1)
    #print var, oadev_x, var/pow(oadev_x,2)
    #(mtaus,mdevs,errs,ns) = at.mdev(x,taus=[af*rate], rate=rate)
    #mdev_x = mdevs[0]
    #rn = pow(mdev_x/oadev_x,2)
    return var/pow(oadev_x,2.0)

def b1(N, mu):
    """
        The exponents are defined as
        S_y(f) = h_a f^alpha
        S_x(f) = g_b f^b
        bias = const * tau^mu
        
        
        and relate to eachother by:
        b    alpha   mu
        0    +2      -2
       -1    +1      -2   resolve between -2 cases with R(n)
       -2     0      -1
       -3    -1       0
       -4    -2      +1
       -5    -3      +2
       -6    -4      +3 for HDEV, by applying B1 to frequency data, and add +2 to resulting mu
    """
    if mu == 2:
        return float(N)*(float(N)+1.0)/6.0
        #up = N*(1.0-pow(N, mu))
        #down = 2*(N-1.0)*(1-pow(2.0, mu))
        #return up/down
        
    elif mu == 1:
        return float(N)/2.0
    elif mu == 0:
        return N*np.log(N)/(2.0*(N-1.0)*np.log(2))
    elif mu == -1:
        return 1
    elif mu == -2:
        return (pow(N,2)-1.0)/(1.5*N*(N-1.0))
    else:
        up = N*(1.0-pow(N, mu))
        down = 2*(N-1.0)*(1-pow(2.0, mu))
        return up/down
        
    assert False # we should never get here

def b_to_mu(b):
    a = b+2
    if a==+2:
        return -2
    elif a==+1:
        return -2
    elif a==0:
        return -1
    elif a==-1:
        return 0
    elif a==-2:
        return 1
    elif a==-3:
        return 2
    elif a==-4:
        return 3
    assert False

def rn_noise_id(x, af, rate):
    """ B1(n) noise identification algorithm
    
    """
    #print "len(x) ", len(x), "oadev",af*rate
    (taus,devs,errs,ns) = at.oadev(x,taus=[af*rate], rate=rate) 
    #print devs
    oadev_x = devs[0]
    (mtaus,mdevs,errs,ns) = at.mdev(x,taus=[af*rate], rate=rate)
    mdev_x = mdevs[0]
    rn = pow(mdev_x/oadev_x,2)
    return rn

ng = at.Noise()
nr=pow(2,14)
qd=2e-20
b=-1

def average_b1(N=20, b=0, af=1):
    ng.set_input(nr,qd,b)
    rns=[]
    for n in range(N):
        ng.generateNoise()
        x = ng.time_series
        #print at.oadev(x,taus=[8000], rate=1) 
        #print af
        rn = b1_noise_id(x, af, 1)
        rns.append(rn)
    return np.average(rns)
    
def average_rn(N=20, b=0, af=1):
    ng.set_input(nr,qd,b)
    rns=[]
    for n in range(N):
        ng.generateNoise()
        x = ng.time_series
        #print at.oadev(x,taus=[8000], rate=1) 
        #print af
        rn = rn_noise_id(x, af, 1)
        rns.append(rn)
    return np.average(rns)

def rn_vs_tau(N=20,b=0):
    rs=[]
    taus=np.logspace(0,np.log( pow(2,5) ),20)
    #print taus
    for t in taus:
        r = average_rn(N,b,int(t))
        rs.append(r)
    return (taus,rs)

def b1_vs_tau(N=20,b=0):
    rs=[]
    taus=np.logspace(0,np.log( pow(2,5) ),20)
    taus = [int(t) for t in taus]
    #print taus
    for t in taus:
        r = average_b1(N,b,int(t))
        rs.append(r)
    return (taus,rs)

def b1_boundary(b_hi, N):
    b_lo = b_hi-1
    b1_lo = b1(N, b_to_mu(b_lo))
    b1_hi = b1(N, b_to_mu(b_hi))
    if b1_lo >= -4:
        return np.sqrt(b1_lo*b1_hi) # geometric mean
    else:
        return 0.5*(b1_lo+b1_hi)

plt.figure()
for b in [0,-1, -2, -3, -4, -5]: # -1,-2,-3,-4
    #(taus, rs) = rn_vs_tau(20,b)
    (btaus, bs) = b1_vs_tau(60,b)
    print b
    plt.loglog(btaus,bs,'o',label="Rn b=%d"%b)
    b1v=[]
    b1_bound=[]
    for t in btaus:
        tmp = range(0,nr,int(t))
        print t, len(tmp)
        b1v.append( b1(len(tmp), b_to_mu(b)) )
        b1_bound.append( b1_boundary(b, len(tmp) ) )
    plt.loglog(btaus,b1v,'-',label="b1() b=%d"%b)
    plt.loglog(btaus,b1_bound,'--',label="boundary(%d, %d)"%(b, b-1))
    
    #plt.loglog(btaus,bs,'o',label="B1 b=%d"%b)

plt.legend()
plt.grid()
plt.show()
