# -*- coding: utf-8 -*-
from pylab import *

### processing function
def store(var,textFile):
    data=loadtxt(textFile,skiprows=1)
    it=[]
    t=[]
    e1=[]
    e2=[]
    e3=[]
    s1=[]
    s2=[]
    s3=[]
    p=[]
    p32=[]
    p33=[]
    tc=[]
    sc=[]
    c=[]
    uf=[]
    for i in range(0,len(data)):
      it.append(float(data[i,3]))
      t.append(float(data[i,11]))
      e1.append(float(data[i,0]))
      e2.append(float(data[i,1]))
      e3.append(float(data[i,2]))
      s1.append(float(-data[i,8]))
      s2.append(float(-data[i,9]))
      s3.append(float(-data[i,10]))
      p.append(float(data[i,4]))
      tc.append(float(data[i,12]))
      sc.append(float(data[i,7]))
      p32.append(float(data[i,5]))
      p33.append(float(data[i,6]))
      uf.append(float(data[i,13]))
    var.append(it)
    var.append(t)
    var.append(e1)
    var.append(e2)
    var.append(e3)
    var.append(s1)
    var.append(s2)
    var.append(s3)
    var.append(p)
    var.append(tc)
    var.append(sc)
    var.append(p32)
    var.append(p33)
    var.append(uf)
   
## data input
dataFile1='111_10k_injection'
a1=[]
store(a1,dataFile1)

### data plot
rcParams.update({'legend.numpoints':1,'font.size': 20,'axes.labelsize':26,'xtick.major.pad':10,'ytick.major.pad':10,'legend.fontsize':20})
lw=2
ms=10
labelSize=15

figure(0,figsize=(14,10))
ax1=subplot(1,1,1)
xlabel('iteration [-]')
ax1.plot(a1[0],[x/1e6 for x in a1[5]],'-r',linewidth=lw)
ax1.plot(a1[0],[x/1e6 for x in a1[6]],'-g',linewidth=lw)
ax1.plot(a1[0],[x/1e6 for x in a1[7]],'-b',linewidth=lw)
ylabel('sigma [MPa]')
legend(('1','2','3'),loc=2)
ax2=ax1.twinx()
ax2.plot(a1[0],[x*1e3 for x in a1[2]],'--r',linewidth=lw)
ax2.plot(a1[0],[x*1e3 for x in a1[3]],'--g',linewidth=lw)
ax2.plot(a1[0],[x*1e3 for x in a1[4]],'--b',linewidth=lw)
ylabel('eps [millisrain]')
#savefig(dataFile1+'_Sig&Eps.pdf',dpi=1000,format='pdf',transparent=False)

figure(1,figsize=(14,10))
plot(a1[0],[x*1e3 for x in a1[2]],'--r',linewidth=lw)
plot(a1[0],[x*1e3 for x in a1[3]],'--g',linewidth=lw)
plot(a1[0],[x*1e3 for x in a1[4]],'--b',linewidth=lw)
ylabel('eps [millisrain]')
legend(('1','2','3'),loc=2)
axis(ymin=0.)
#savefig(dataFile1+'_Sig&Eps.pdf',dpi=1000,format='pdf',transparent=False)


figure(2,figsize=(14,10))
ax1=subplot(1,1,1)
#axis(xmin=0.,xmax=2.)
xlabel('iteration [-]')
ax1.plot(a1[0],[x/1e6 for x in a1[8]],'-r',linewidth=1.5*lw)
ylabel('p [MPa]',color='r')
ax2=ax1.twinx()
#axis(xmin=0.,xmax=2.)
ax2.plot(a1[0],[i+j for i,j in zip(a1[9],a1[10])],'b',linewidth=lw)
ylabel('microcracks [-]',color='b')
#axis(ymin=0.,ymax=1500)
#frame1 = plt.gca()
#frame1.axes.yaxis.set_ticklabels([])
#frame1.axes.yaxis.set_ticks(arange(0,1500.01,300))
savefig(dataFile1+'_pressure&Cracks.pdf',dpi=1000,format='pdf',transparent=False)

#figure(3,figsize=(8,6))
#grid()
#ax1=subplot(1,1,1)
#xlabel('iteration [-]')
#ax1.plot(a1[0],[x/1e6 for x in a1[8]],'-r',linewidth=1.5*lw)
#ylabel('p [psi]',color='r')
#ax2=ax1.twinx()
#ax2.plot(a1[0],a1[12],'b',linewidth=0.5*lw)
#ylabel('unb force [-]',color='b')
#axis(ymin=0.,ymax=0.1)
#savefig(dataFile1+'_pressure&unbF.pdf',dpi=1000,format='pdf',transparent=False)

### show or save
show()

