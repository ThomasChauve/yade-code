# -*- coding: utf-8 -*-
from pylab import *

### processing function
def store(var,textFile):
    data=loadtxt(textFile,skiprows=1)
    it=[]
    t=[]
    e=[]
    e1=[]
    e2=[]
    s=[]
    tc=[]
    sc=[]
    unbF=[]
    for i in range(5,len(data)):
      it.append(float(data[i,3]))
      e.append(float(data[i,0]))
      e1.append(float(data[i,2]))
      e2.append(float(data[i,1]))
      s.append(float(data[i,5]))
      tc.append(float(data[i,6]))
      sc.append(float(data[i,4]))
      unbF.append(float(data[i,7]))
    var.append(it)
    var.append(e)
    var.append(e1)
    var.append(e2)
    var.append(s)
    var.append(tc)
    var.append(sc)
    var.append(unbF)

### data input
dataFile='010201_5k_3_K10_tensionTest_r0.02'
a=[]
store(a,dataFile)

### data plot
rcParams.update({'font.size':28,'axes.labelsize':36,'axes.labelpad':5,'xtick.major.pad':10,'ytick.major.pad':10,'legend.fontsize':28,'legend.numpoints':1})

lw=2
ms=14

figure(1,figsize=(11,11))
ax1=subplot(1,1,1)
xlabel('iteration [-]')
ax1.plot(a[0],[i/1e6 for i in a[4]],'-k',linewidth=lw)
ylabel(r'$\sigma_1$ [MPa]')
#legend((r'$\sigma_1$'),loc=3)
ax2=ax1.twinx()
ax2.plot(a[0],[i*1e3 for i in a[1]],'--r',linewidth=lw)
ax2.plot(a[0],[i*1e3 for i in a[2]],'-r',linewidth=lw)
ax2.plot(a[0],[i*1e3 for i in a[3]],'-g',linewidth=lw)
ylabel(r'$\varepsilon_i$ [millistrain]')
legend((r'$\varepsilon_1$ ref',r'$\varepsilon_1$',r'$\varepsilon_2$'),loc=2)
#savefig('SigVsIter.tiff',dpi=200,format='tiff',transparent=False)

figure(2,figsize=(14,11))
ax1=subplot(1,1,1)
grid()
xlabel(r'$\varepsilon_1$ [millistrain]')
ax1.plot([abs(i)*1e3 for i in a[1]],[i/j/1e9 for i,j in zip(a[4],a[1])],'-r',linewidth=lw)
ylabel(r'$E$ [GPa]',color='red')
#axis(xmax=20000,ymin=18,ymax=26)
axis(xmax=0.1,ymin=21,ymax=21.5)
ax2=ax1.twinx()
grid()
ax2.plot([abs(i)*1e3 for i in a[1]],[-i/j for i,j in zip(a[3],a[1])],'-b',linewidth=lw)
ylabel(r'$\nu$ [-]',color='blue')
#axis(xmax=20000,ymin=0.1,ymax=0.3)
axis(xmax=0.1,ymin=0.13,ymax=0.17)
#savefig(dataFile+'_E&NuVsIter.tiff',dpi=200,format='tiff',transparent=False)

figure(3,figsize=(14,11)) #,figsize=(14,10)
ax1=subplot(1,1,1)
xlabel(r'$\varepsilon_1$ [millistrain]')
ax1.plot([abs(i)*1e3 for i in a[1]],[abs(i)/1e6 for i in a[4]],'-b',linewidth=lw)
ylabel(r'$\sigma_1\ \mathrm{[MPa]}$',color='blue')
axis(xmax=0.4)
ax2=ax1.twinx()
ax2.plot([abs(i)*1e3 for i in a[1]],[(i+j) for i, j in zip(a[5],a[6])],'-r',linewidth=lw)
ylabel('$number\ of\ microcracks\ \mathrm{[-]}$',color='red')
axis(xmax=0.4)
ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
#savefig(dataFile+'_Sig&CracksVsEps.tiff',dpi=500,format='tiff',transparent=False)
#savefig(dataFile+'_Sig&CracksVsEps.ps',dpi=400,format='ps',transparent=False)

print 'UTS=',max(a[4])/1.e6

show()
