"""

##########################
#  test stable32plot.py  #
##########################

""" 
import allantools
from pylab import figure,show
from stable32plot import sigmaplot,dataplot#import 2 functions: sigmaplot,dataplot


"""#------------generate random data and cal adev-----------------"""
x1 = allantools.noise.white(1000)
(taus, adevs, errors, ns) = allantools.adev(x1)
(taust, adevst, errorst, nst) = allantools.tdev(x1)
adevs2=[]
adevs3=[]
adevs4=[]
adevs5=[]
adevs6=[]
adevsn=[]
errors2=[]
errors3=[]
errors4=[]
errors5=[]
errors6=[]

for i in range(len(taus)):
    adevs2.append(0.6*adevs[i])
    adevs3.append(0.5*adevs[i])
    adevs4.append(0.4*adevs[i])
    adevs5.append(0.3*adevs[i])
    adevs6.append(0.2*adevs[i])
    adevsn.append(0.05*adevs[i])
    errors2.append(0.6*errors[i])
    errors3.append(0.5*errors[i])
    errors4.append(0.4*errors[i])
    errors5.append(0.3*errors[i])
    errors6.append(0.2*errors[i])
    
x2 = allantools.noise.white(500)

"""#--------------eg1: plot single adev curve with text list----------------"""
figure()
sigmaplot(taus,adevs,errors)

"""#--------------eg2: plot single adev curve with noise floor of the test device--------"""
figure()
sigmaplot(taus,adevs,errors,"test 1",taun=taus,sigman=adevsn)

"""#--------------eg3: plot single tdev curve with text list----------------"""
figure()
sigmaplot(taust,adevst,errorst,sigmatype="tdev")

"""#--------------eg4: plot multiple adev curve with noise floor of the test device"""
"""# eg4 contains all the args """ 
figure()
sigmaplot(taus,adevs,errors,"test 1",
         tau2=taus,sigma2=adevs2,error2=errors2,legend2="test 2",
         tau3=taus,sigma3=adevs3,error3=errors3,legend3="test 3",
         tau4=taus,sigma4=adevs4,error4=errors4,legend4="test 4",
         tau5=taus,sigma5=adevs5,error5=errors5,legend5="test 5",
         tau6=taus,sigma6=adevs6,error6=errors6,legend6="test 6",
         taun=taus,sigman=adevsn,legendn="test 1 noise floor",
         sigmatype="oadev",sigmatext=False)

"""#-----------------eg5: plot single frequency data-----------"""
figure()
dataplot(range(len(x2)),x2)

"""#--------------eg6: plot multiple phase data----------------"""
figure()
dataplot(range(len(x1)),x1,"test 1",
         sec2=range(len(x2)),data2=x2,legend2="test 2",
         datatype='phase')

show()


