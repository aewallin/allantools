"""
############################################
#    stable32 style sigma plot function    #
#                                          #
#   Author: Yan Xie (y.xie_cn@outlook.com) #
############################################

1)function sigmaplot:
--generate stable32 style *DEV plot;
--plot maximum 6 curves with noise floor;
--display tau-sigma list, setting 'sigmatext=False' to hide the text list;
--supported sigmatype includes:
  adev,oadev,madev,tdev,hdev,ohdev,totdev,ttotdev,mtotdev,htotdev
-- **arg keywords include:
    tau2,sigma2,error2,legend2: plot the 2nd curve
    tau3,sigma3,error3,legend3: plot the 3rd curve
    tau4,sigma4,error4,legend4: plot the 4th curve
    tau5,sigma5,error5,legend5: plot the 5th curve
    tau6,sigma6,error6,legend6: plot the 6th curve
    taun,sigman,errorn: plot noise floor
    sigmatype: supported sigmaplot types includes 'adev','oadev','madev','tdev','hdev','ohdev',
               'totdev','ttotdev','mtotdev','htotdev' (default: 'adev')
    sigmatext: is a bool variable, show/hide tau1-sigma1 text list in the plot
                (default: True)

2) function dataplot:
--plot stable32 style frequency/phase data;
--plot maximum 6 curves;
--**arg keywords include:
    sec2,data2,legend2: plot the 2nd curve
    sec3,data3,legend3: plot the 3rd curve
    sec4,data4,legend4: plot the 4th curve
    sec5,data5,legend5: plot the 5th curve
    sec6,data6,legend6: plot the 6th curve
    datatype: supported data types includes 'freq' and 'phase'(default: 'freq')

""" 


from pylab import *

"""#--------------function plot sigma ---------begin-----------"""
def sigmaplot(tau1,sigma1,error1,legend1="",**arg):
    
    """#all the supported sigmatype tuple"""
    dev_tuple=('adev','oadev','mdev','tdev','hdev','ohdev',
              'totdev','ttotdev','mtotdev','htotdev')
    
    """#setting plot sigmatype"""
    sigmatype="adev"
    if 'sigmatype' in arg.keys():
        if arg['sigmatype'] in dev_tuple:
            sigmatype=arg['sigmatype']
        #else:
            #add an "error" report
        
    """#setting show/hide tau-sigma list text"""
    sigmatext=True    
    if 'sigmatext'in arg.keys():
        sigmatext=arg['sigmatext']

    """#set curve color tuple"""
    color1='#B30000' #dark red
    color2='#33A1C9'#blue
    color3='#CD2990'#maroon
    color4='#32CD32'#green
    color5='#FF8000'#orange
    color6='#8E388E'#purple
    colorn='#C9C9C9' #gray
    
    """#set arg key tuple"""
    taui_tuple=('tau2','tau3','tau4','tau5','tau6')
    sigmai_tuple=('sigma2','sigma3','sigma4','sigma5','sigma6')
    errori_tuple=('error2','error3','error4','error5','error6')
    legendi_tuple=('legend2','legend3','legend4','legend5','legend6')
    colori_tuple=(color2,color3,color4,color5,color6)
    
    """#plot errorbars"""
    ax=axes()
    ax.set_xscale("log")
    ax.set_yscale("log")
    fillstyle1=set_markerfillstyle(sigmatype)
    """#plot 1st curve"""
    leghandle1=ax.errorbar(tau1,sigma1,yerr=error1,capsize=3,color=color1,
                marker='.',fillstyle=fillstyle1,markersize=12)
    
    """#set variable list for the ith curve"""
    taui=[]
    sigmai=[]
    errori=[]
    leghandlei=[]
    legendi=[]

    """#plot ith curve"""
    for i in range(len(taui_tuple)):
        if (taui_tuple[i] in arg.keys()) and (sigmai_tuple[i] in arg.keys()) and (errori_tuple[i] in arg.keys()):
            taui.append(arg[taui_tuple[i]])
            sigmai.append(arg[sigmai_tuple[i]])
            errori.append(arg[errori_tuple[i]])
            leghandlei.append(ax.errorbar(taui[len(taui)-1],sigmai[len(sigmai)-1],
                                    yerr=errori[len(errori)-1],capsize=3,
                                    color=colori_tuple[i],marker='.',
                                    fillstyle=fillstyle1,markersize=12))
            if legendi_tuple[i] in arg.keys():
                legendi.append(arg[legendi_tuple[i]])
            else:
                legendi.append("")

    """#plot noise floor"""
    if ('taun' in arg.keys()) and ('sigman' in arg.keys()):
            leghandlen=ax.fill_between(arg['taun'],arg['sigman'],color=colorn)
            if 'legendn' in arg.keys():
                legendn=arg['legendn']
            else:
                legendn="noise floor"
    else:
        legendn=""

    """#show legend"""
    if sigmatext==True: #set legend location
        loci=3
    else:
        loci=0
    if legend1!="":
        legendi.insert(0,legend1)
        leghandlei.insert(0,leghandle1)
        if legendn!="":
            legendi.insert(len(legendi),legendn)
            leghandlei.insert(len(leghandlei),leghandlen)
        ax.legend(leghandlei,legendi,loc=loci)
                   
    """#show grid&ticks"""
    show_grid()
    
    """#show tau&sigma list text"""
    xlimval=ax.get_xlim()
    ylimval=ax.get_ylim()
    tlocx=xlimval[1]-0.12*(xlimval[1]-xlimval[0])#set text x-axis location
    tlocy=ylimval[1]-0.12*(ylimval[1]-ylimval[0])#set text y-axis location
    if sigmatext==True:
        if len(tau1)<=20:
            show_text(tlocx,tlocy,tau1,sigma1)
        else:# to be changed
            disptau1=tau1[0:20]
            dispsigma1=sigma1[0:20]
            show_text(tlocx,tlocy,disptau1,dispsigma1)
                      
    """#show xlabel,ylabel,title"""
    show_label_title(sigmatype)
"""#--------------function plot sigma ---------end-----------"""


def set_markerfillstyle(tmpsigmatype):    
    if tmpsigmatype=='tdev' or tmpsigmatype=='ttotdev': 
        tmpfillstyle='none'
    else:
        tmpfillstyle='full'
    return tmpfillstyle
        

def show_grid():
    minorticks_on()
    grid(which='major',linestyle='-')
    grid(which='minor',linestyle='--')
    xticks(fontsize=12,fontweight='bold',family='Times New Roman')
    yticks(fontsize=12,fontweight='bold',family='Times New Roman',rotation=90)

def show_text(lx,ly,tmptau,tmpsigma):
    disp_text="   Tau"+'{0:11}'.format(' ')+"Sigma \n"
    for i in range(len(tmptau)):
        disp_text+='{0:7.2e}'.format(tmptau[i])+'{0:4}'.format(' ')+'{0:7.2e}'.format(tmpsigma[i])
        if i<=(len(tmptau)-2):
            disp_text+='\n'   
    text(lx,ly,disp_text,size=9,
        va="top", ha="right", multialignment="center",
        family='Century Gothic',
        bbox=dict(fc='white',linewidth=0.5))

def show_label_title(tmpsigmatype):
    """#set ylabel"""
    if tmpsigmatype=="adev":
        str_ylabel="Allan Deviation, \u03C3"+'$_y$'+"(\u03C4)"
    elif tmpsigmatype=="tdev":
        str_ylabel="Time Deviation, \u03C3"+'$_x$'+"(\u03C4), Seconds"
    elif tmpsigmatype=="oadev":
        str_ylabel="Overlapping Allan Deviation, \u03C3"+'$_y$'+"(\u03C4)"
    elif tmpsigmatype=="mdev":
        str_ylabel="Modified Allan Deviation, Mod \u03C3"+'$_y$'+"(\u03C4)"
    elif tmpsigmatype=="hdev":
        str_ylabel="Hadamard Deviation, H\u03C3"+'$_y$'+"(\u03C4)"
    elif tmpsigmatype=="ohdev":
        str_ylabel="Overlapping Hadamard Deviation, H\u03C3"+'$_y$'+"(\u03C4)"
    elif tmpsigmatype=="totdev":
        str_ylabel="Total Allan Deviation, \u03C3"+'$_{total}$'+"(\u03C4)"
    elif tmpsigmatype=="ttotdev":
        str_ylabel="Total Time Deviation, \u03C3"+'$_x$'+"(\u03C4), Sec"
    elif tmpsigmatype=="mtotdev":
        str_ylabel="Total Modified Deviation, Mod \u03C3"+'$_{total}$'+"(\u03C4)"
    elif tmpsigmatype=="htotdev":
        str_ylabel="Total Hadamard Deviation, H\u03C3"+'$_{total}$'+"(\u03C4)"
    else:
        str_ylabel="Allan Deviation, \u03C3"+'$_y$'+"(\u03C4)"
        
    """#set title """   
    if tmpsigmatype=='tdev' or tmpsigmatype=='ttotdev':
        str_title="TIME STABILITY"
    else:
        str_title="FREQUENCY STABILITY"
    
    """#display xlabel,ylabel,title"""
    xlabel("Averaging Time, \u03C4, Seconds",fontsize=16,
           color='#00008B',family='Times New Roman',fontweight='bold')
    ylabel(str_ylabel,
           family='Times New Roman',fontsize=16,color='#00008B',fontweight='bold')
    title(str_title,fontsize=18,
          family='Arial',color='#00008B',fontweight='bold')

    
"""#-------function plot freq/phase data-------------begin---------------"""
def dataplot(sec1,data1,legend1="",**arg):
    """#all the supported data type tuple"""
    datatype_tuple=('freq','phase')
    
    """#setting plot datatype"""
    datatype="freq"
    if 'datatype' in arg.keys():
        if arg['datatype'] in datatype_tuple:
            datatype=arg['datatype']
            
    """#set curve color tuple"""
    color1='#B30000' #dark red
    color2='#33A1C9'#blue
    color3='#CD2990'#maroon
    color4='#32CD32'#green
    color5='#FF8000'#orange
    color6='#8E388E'#purple
    """#set arg key tuple"""
    seci_tuple=('sec2','sec3','sec4','sec5','sec6')
    datai_tuple=('data2','data3','data4','data5','data6')
    legendi_tuple=('legend2','legend3','legend4','legend5','legend6')
    colori_tuple=(color2,color3,color4,color5,color6)


    """#plot curves"""
    """#plot 1st curve"""
    leghandle1,=plot(sec1,data1,color=color1)
    """#set variable list for the ith curve"""
    seci=[]
    datai=[]
    leghandlei=[]
    legendi=[]
    """#plot ith curve"""
    for i in range(len(seci_tuple)):
        if (seci_tuple[i] in arg.keys()) and (datai_tuple[i] in arg.keys()):
            seci.append(arg[seci_tuple[i]])
            datai.append(arg[datai_tuple[i]])
            tmphandle,=plot(seci[len(seci)-1],datai[len(datai)-1],                                   
                             color=colori_tuple[i])
            leghandlei.append(tmphandle)
            if legendi_tuple[i] in arg.keys():
                legendi.append(arg[legendi_tuple[i]])
            else:
                legendi.append("")
    
    """#show legend"""
    if legend1!="":
        legendi.insert(0,legend1)
        leghandlei.insert(0,leghandle1)
        legend(leghandlei,legendi)

    """#show grid & ticks """   
    minorticks_on()
    grid(which='major',linestyle='-')
    grid(which='minor',linestyle=':')
    xticks(fontsize=12,family='Times New Roman',fontweight='bold')
    yticks(fontsize=12,family='Times New Roman',fontweight='bold')
    
    """#display label & title"""
    xlabel("Data Point",fontsize=16,color='#00008B',
           family='Times New Roman',fontweight='bold')
    if datatype=='freq':
        str_ylabel="Frequency"
        str_title="FREQUENCY DATA"
    elif datatype=='phase':
        str_ylabel="Phase, Seconds"
        str_title="PHASE DATA"
    else:#change to an "error" report
        str_ylabel="Frequency"
        str_title="FREQUENCY DATA"   
    ylabel(str_ylabel,fontsize=16,color='#00008B',
           family='Times New Roman',fontweight='bold')
    title(str_title,fontsize=18,color='#00008B',
          family='Arial',fontweight='bold')

"""#-------function plot freq/phase data-------------end--------------- """
        
