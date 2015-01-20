import numpy as np
import matplotlib 
matplotlib.use('PDF')
import pylab as py

sol = 3e5
h=0.67

Z=0.05
print "Z="+str(Z)
if (Z==5):
    class_d = np.genfromtxt('./Testing/class/th2_mnu0p05_z6_tk.dat')
    class_v = np.genfromtxt('./Testing/class/th2_v_mnu0p05_z3_tk.dat')
    class_l = np.genfromtxt('./Testing/class/th2_mnu0p05_z6_v_tk.dat')
    camb_d = np.genfromtxt('./Testing/camb/interp_mnu0p05_z5.dat')
    camb_v = np.genfromtxt('./Testing/camb/vinterp_mnu0p05_z5.dat')
if (Z==10):
    class_d = np.genfromtxt('./Testing/class/th2_mnu0p05_z3_tk.dat')
    class_v = np.genfromtxt('./Testing/class/th2_v_mnu0p05_z3_tk.dat')
    class_l = np.genfromtxt('./Testing/class/th2_mnu0p05_z3_v_tk.dat')
    camb_d = np.genfromtxt('./Testing/camb/sim_mnu0p05_transfer_out_z10.dat')
    camb_v = np.genfromtxt('./Testing/camb/sim_mnu0p05_veltransfer_out_z10.dat')
if (Z==0):
    class_d = np.genfromtxt('./Testing/class/th2_mnu0p05_z13_tk.dat')
    class_v = np.genfromtxt('./Testing/class/th2_v_mnu0p05_z13_tk.dat')
    class_l = np.genfromtxt('./Testing/class/th2_mnu0p05_z12_v_tk.dat')
    camb_d = np.genfromtxt('./Testing/camb/sim_mnu0p05_transfer_out_z0.dat')
    camb_v = np.genfromtxt('./Testing/camb/sim_mnu0p05_veltransfer_out_z0.dat')
if (Z==0.05):
    class_d = np.genfromtxt('./Testing/class/th2_mnu0p05_z12_tk.dat')
    class_v = np.genfromtxt('./Testing/class/th2_v_mnu0p05_z12_tk.dat')
    class_l = np.genfromtxt('./Testing/class/th2_mnu0p05_z12_v_tk.dat')
    camb_d = np.genfromtxt('./Testing/camb/sim_mnu0p05_transfer_out_z0p05.dat')
    camb_v = np.genfromtxt('./Testing/camb/sim_mnu0p05_veltransfer_out_z0.dat')
    
#mink = class_d[0,0]
#maxk = class_d[class_d.shape[0]-1,0]
mink = camb_d[0,0]
maxk = 200 #camb_d[camb_d.shape[0]-1,0]
Nk = 1000

eps = 1e-7
newk = np.logspace(np.log10(mink+eps),np.log10(maxk-eps),Nk)

def lin_interp(tf,k):
    rtf = np.zeros((k.size,tf.shape[1]))
    rtf[:,0]=k
    for i in range(tf.shape[1]-1):
        rtf[:,i+1] = np.interp(k,tf[:,0],tf[:,i+1])
    return rtf

def log_interp(tf,k):
    rtf = np.zeros((k.size,tf.shape[1]))
    rtf[:,0]=k
    for i in range(tf.shape[1]-1):
        rtf[:,i+1] = 10.0**np.interp(np.log10(k),np.log10(tf[:,0]),np.log10(tf[:,i+1]))
    return rtf

class_di = log_interp(class_d,newk)
class_vi = log_interp(np.abs(class_v),newk)
class_vi[:,1:class_vi.shape[0]-1] = -1.0*class_vi[:,1:class_vi.shape[0]-1]
class_li = log_interp(np.abs(class_l),newk)
class_li[:,1:class_vi.shape[0]-1] = -1.0*class_li[:,1:class_vi.shape[0]-1]

#class_dc = log_interp(class_d,camb_d[:,0])
#class_vc = log_interp(class_v,camb_v[:,0])

#Write out interpolated files
if (Z==10):
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z10_tk.dat',class_di)
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z10_v_tk.dat',class_li)
if (Z==5):
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z5_tk.dat',class_di)
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z5_v_tk.dat',class_li)
if (Z==0):
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z0_tk.dat',class_di)
    np.savetxt('./Testing/Interp/ith2_v_mnu0p05_z0_tk.dat',class_vi)
if (Z==0.05):
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z0p05_tk.dat',class_di)
    np.savetxt('./Testing/Interp/ith2_mnu0p05_z0p05_v_tk.dat',class_li)
#Matplotlib plot
py.figure(1)

py.loglog(camb_d[:,0],camb_d[:,1],'--',label='CAMB DM',color='blue')
py.loglog(camb_d[:,0],camb_d[:,5],'--',label='CAMB NU',color='red')

py.loglog(class_d[:,0],class_d[:,1],':',label='CLASS DM',color='blue')
py.loglog(class_d[:,0],class_d[:,5],':',label='CLASS NU',color='red')

py.loglog(class_di[:,0],class_di[:,1],'-',label='CLASS DM INTERP',color='blue')
py.loglog(class_di[:,0],class_di[:,5],'-',label='CLASS NU INTERP',color='red')

xmin = mink
xmax = maxk

py.xlabel('k (h/Mpc)')
py.axes().set_xlim((xmin,xmax))
py.ylabel(r'TF(k)')
#py.axes().set_ylim((ymin,ymax))
py.legend(loc='best',fancybox=True)
py.savefig('./Testing/Interp_d.pdf')
py.close()

#Matplotlib plot
py.figure(1)

py.loglog(camb_v[:,0],np.abs(camb_v[:,1]),':',label='CAMB DM',color='blue')
py.loglog(camb_v[:,0],np.abs(camb_v[:,5]),':',label='CAMB NU',color='red')

py.loglog(class_v[:,0],sol*np.abs(class_v[:,1])/(class_v[:,0]*h)**3,':',label='CLASS DM',color='grey')
py.loglog(class_v[:,0],sol*np.abs(class_v[:,5])/(class_v[:,0]*h)**3,':',label='CLASS NU',color='grey')

py.loglog(class_vi[:,0],sol*np.abs(class_vi[:,1])/(class_vi[:,0]*h)**3,'--',label='CLASS DM INTERP',color='grey')
py.loglog(class_vi[:,0],sol*np.abs(class_vi[:,5])/(class_vi[:,0]*h)**3,'--',label='CLASS NU INTERP',color='grey')

py.loglog(class_li[:,0],np.abs(class_li[:,1]),'--',label='CLASS DM LINEAR INTERP',color='blue')
py.loglog(class_li[:,0],np.abs(class_li[:,5]),'--',label='CLASS NU LINEAR INTERP',color='red')

print "Vel ratio class lin to class newt DM = " + str(class_li[0,1]/(sol*np.abs(class_vi[0,1])/(class_vi[0,0]*h)**3))
print "Vel ratio class lin to class newt NU = " + str(class_li[0,5]/(sol*np.abs(class_vi[0,5])/(class_vi[0,0]*h)**3))

camb_vi = log_interp(np.abs(camb_v),newk)
camb_vi[:,1:camb_vi.shape[0]-1] = -1.0*camb_vi[:,1:camb_vi.shape[0]-1]

print "Vel ratio class lin to camb lin DM = " + str(class_li[0,1]/camb_vi[0,1])
print "Vel ratio class lin to camb lin NU = " + str(class_li[0,5]/camb_vi[0,5])

xmin = mink
xmax = maxk

py.xlabel('k (h/Mpc)')
py.axes().set_xlim((xmin,xmax))
py.ylabel(r'TF(k)')
#py.axes().set_ylim((ymin,ymax))
py.legend(loc='best',fancybox=True)
py.savefig('./Testing/Interp_v.pdf')
py.close()


