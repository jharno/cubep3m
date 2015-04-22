#!/usr/bin/python
import numpy as np

doPlot = False

#z of fE
Z = 10
ZE = Z
if (ZE == 10):
    ZA = 10.1
    ZB = 9.9
    ZE = 10
    print "Computing at " + str(ZE)


if (ZE == 0):
    ZA = 0.1
    ZB = 0.0
    ZE = 0.05
    print "Computing at " + str(ZE) + " as proxy for " + str(0)
#Location of transfer functions

mnu='3mnu0p2'

#dir = '/Users/derekinman/Research/Camb/mnu1p0/linear/'
#dir = '/Users/derekinman/Research/Camb/mnu0p2/non_linear/'
#prefix = 'nu_'+mnu+'_transfer_out_z'
#suffix = '.dat'

dir = './3mnu0p2/' #/Users/derekinman/Research/Camb/simulation/mnu0p05/'
prefix = 'sim_3mnu0p2_transfer_out_z'
suffix = '.dat'

if (Z == 10):
    fE = prefix+str(ZE)+suffix
    fA = prefix+str(int(ZA))+'p1'+suffix
    fB = prefix+str(int(ZB))+'p9'+suffix
    fO = 'sim_'+mnu+'_veltransfer_out_z'+str(ZE)+suffix
if (Z == 0):
    fE = prefix+'0p05'+suffix
    fA = prefix+'0p1'+suffix
    fB = prefix+'0'+suffix
    fO = 'sim_3mnu0p2_veltransfer_out_z0'+suffix
#scale factor
a = 1.0/(1.0+ZE)

#Cosmological parameters
OM = 0.32
OL = 1.0 - OM

#Assuming k-col in h/Mpc
H0 = 100

#Compute H(a)
h = 0.67
H = H0 * ( OM*a**-3 + OL )**0.5

#Get transfer functions into arrays
tfE = np.genfromtxt(dir+fE)
tfA = np.genfromtxt(dir+fA)
tfB = np.genfromtxt(dir+fB)

if(tfA.shape!=tfB.shape or tfA.shape!=tfE.shape):
    print "ERROR - TRANSFER FUNCTIONS DIFFERENT SHAPE"
    
Nk = tfA.shape[0]
Nc = tfA.shape[1]

#Check if transfer functions have same k-values
if( not( np.array_equal( tfA[:,0], tfB[:,0] ) and np.array_equal( tfA[:,0], tfE[:,0] ) ) ):
    print "ERROR - TRANSFER FUNCTIONS HAVE DIFFERENT K-VALUES"
    
tfV = np.zeros((Nk,Nc))
for k in range(Nk):
    for c in range(Nc):
        if (c==0):
            tfV[k,c] = tfE[k,c]
        else:
            tfV[k,c] = (H/tfE[k,0])*(tfA[k,c]-tfB[k,c])/(ZA-ZB) #proper km/s
    
np.savetxt(dir+fO,tfV)

#Flatness plot:
if doPlot:
    F = (tfV[:,1]/tfE[:,1])/(-1.0*a*H*(OM*a**-3/(OM*a**-3+OL))**0.6/tfE[:,0])   
    Fplot = list_plot( zip(tfV[:,0],F), plotjoined=true )
    Fplot.show()

    #Test plot
    for k in range(Nk):
        for c in range(Nc):
            if (c==0):
                tfV[k,c] = tfE[k,c]
            else:
                tfV[k,c] = tfV[k,c] * tfE[k,0]**2 * h**2 
            
    D=2.42*10**-9
    spd = 1#299792.458
    plot = list_plot(zip(tfV[:,0],tfV[:,1]**2*D/spd**2),plotjoined=true,color='blue')
    plot+= list_plot(zip(tfV[:,0],tfV[:,2]**2*D/spd**2),plotjoined=true,color='red')
    plot+= list_plot(zip(tfV[:,0],(tfV[:,2]-tfV[:,1])**2*D/spd**2),plotjoined=true,color='orange')
    plot.show(scale='loglog')
    
    plot = list_plot(zip(tfV[:,0],(tfV[:,2]-tfV[:,1])**2*D/spd**2),plotjoined=true,color='orange')
    plot.show(scale='semilogx')
    
    
