# import modules
from numpy import *
from matplotlib import pyplot as plt

import os # directory fucntions
import csv # for reading data
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit # fitting algorithms



# TWO LAYER MODEL

#            air
#-------------------------
#       N1 (thin film)
#-------------------------
#       N2 (substrate)

# Fresnel equations for first interface

def r01_p(N1, theta_0):

    Nti=N1 # = Nt/Ni
    a= Nti**2 * cos(theta_0) - sqrt(Nti**2 - (sin(theta_0))**2)
    b= Nti**2 * cos(theta_0) + sqrt(Nti**2 - (sin(theta_0))**2)
    return a/b

def r01_s(N1, theta_0):

    Nti=N1 # = Nt/Ni
    a= cos(theta_0) - sqrt(Nti**2 - (sin(theta_0))**2)
    b= cos(theta_0) + sqrt(Nti**2 - (sin(theta_0))**2)
    return a/b

# test plot, rp should go to 0 at Brester angle when N are real numbers, air/glass (n=1.33) interface
n_glass=1.46
angles = linspace(deg2rad(0),deg2rad(90),500)
data_rp = absolute(r01_p(n_glass + 0j, angles)) # interesting: observe how adding absorption (k>0) increases reflectivity at Brewster angle
data_rs = absolute(r01_s(n_glass + 0j, angles))
    
plt.cla()
plt.axis([0,90,0,1])
plt.plot(rad2deg(angles), data_rp,'b-', label='|rp_1|')
plt.plot(rad2deg(angles), data_rs,'r-', label='|rs_1|')


plt.axhline(0, color='black')
plt.axvline(0, color='black')

plt.legend()
plt.show()


# ## Fresnel equations for second interface


def r12_p(N1, N2, theta_0):
    theta_1 = arccos((1/N1) * sqrt(1 - (1/N1)**2 * (sin(theta_0))**2)) # first interface transmission angle
    #cos_theta_2 = sqrt(1 - (N1/N2)**2 * (sin(theta_1))**2)
    #theta_2 = arccos(cos_theta_2)
    
    Nti=N2/N1 # = Nt/Ni
    a= Nti**2 * cos(theta_1) - sqrt(Nti**2 - (sin(theta_1))**2)
    b= Nti**2 * cos(theta_1) + sqrt(Nti**2 - (sin(theta_1))**2)
    return a/b

def r12_s(N1, N2, theta_0):
    theta_1 = arccos((1/N1) * sqrt(1 - (1/N1)**2 * (sin(theta_0))**2)) # first interface transmission angle

    Nti=N2/N1 # = Nt/Ni
    a= cos(theta_1) - sqrt(Nti**2 - (sin(theta_1))**2)
    b= cos(theta_1) + sqrt(Nti**2 - (sin(theta_1))**2)
    return a/b

# test plot
plt.cla()

angles=linspace(deg2rad(0),deg2rad(90),500)

N1_test = 1.0 + 0.0j
N2_test = 1.46 + 0.0j # test with air/glass interface
plt.axis([0,90,0,1])
plt.plot(rad2deg(angles),absolute(r12_p(N1_test,N2_test, angles)),'b-', label='|rp_2|') # we plot absolute value of COMPLEX reflection amplitude
plt.plot(rad2deg(angles),absolute(r12_s(N1_test,N2_test, angles)),'r-', label='|rs_2|')

plt.legend()
plt.show()


# ## Fresnel equations for all layers

lambd= 658

def r012_p(lamd, d, N1, N2, theta_0):
    beta = ((2*pi*d)/lambd) * sqrt(N1**2 - 1**2 * (sin(theta_0))**2)
    num = r01_p(N1, theta_0) + r12_p(N1,N2, theta_0) * exp(-1j*2*beta)
    denum = 1 + r01_p(N1, theta_0) * r12_p(N1, N2, theta_0) * exp(-1j*2*beta)
    return num/denum

def r012_s(lamd, d, N1, N2, theta_0):
    beta = ((2 * pi * d)/lambd) * sqrt(N1**2 - 1**2 * (sin(theta_0))**2)
    num = r01_s(N1, theta_0) + r12_s(N1, N2, theta_0) * exp(-1j*2*beta)
    denum = 1 + r01_s(N1, theta_0) * r12_s(N1, N2, theta_0) * exp(-1j*2*beta)
    return num/denum

# test plot
plt.cla()

angles=linspace(deg2rad(0),deg2rad(90),500)

plt.axis([0,90,0,1])
plt.plot(rad2deg(angles),absolute(r012_p(lambd, 10, 1.33, 2.1, angles)),'b-', label='|r012_p|')
plt.plot(rad2deg(angles),absolute(r012_s(lambd, 10, 1.33, 2.1, angles)),'r-', label='|r012_s|')
plt.legend(loc='lower left')
plt.show()


## ellipsometric quantities
def psi(lamd, d, N1, N2, theta_0):
    rp_abs = absolute(r012_p(lamd, d, N1, N2, theta_0))
    rs_abs = absolute(r012_s(lamd, d, N1, N2, theta_0))
    return arctan2(rp_abs, rs_abs) #see the definition of 'arctan2', it's reverse of 'ordinary' -> arctan2(y,x)

def delta(lamd, d, N1, N2, theta_0):
    rp = r012_p(lamd, d, N1, N2, theta_0)
    rs = r012_s(lamd, d, N1, N2, theta_0)
    return angle(rp/rs)


# Test plots
# test plot psi

#     air
#---------------
#    N1: SiO2
#---------------
#  N2: Si wafer

plt.cla()

angles=linspace(deg2rad(40),deg2rad(80),500)

#plt.axis([0,90,-1,1])
N1plot = 1.48 + 0.0j     #SiO2
N2plot = 3.5 + 0.014j      #Si

plt.cla()
plt.plot(rad2deg(angles), rad2deg(psi(lambd, 0.3, N1plot, N2plot, angles)),'r-', label=r'd(SiO$_2$) = 0.3 nm')
plt.plot(rad2deg(angles), rad2deg(psi(lambd, 5., N1plot, N2plot, angles)),'g-', label=r'd(SiO$_2$) = 5.0 nm')
plt.plot(rad2deg(angles), rad2deg(psi(lambd, 10., N1plot, N2plot, angles)),'b-', label=r'd(SiO$_2$) = 10.0 nm')
plt.xlabel(r'$\theta_{i}$ $(^{\circ})$',size=20)
plt.ylabel(r'$\psi$ $(^{\circ})$',size=25,rotation=90)
plt.legend(loc='lower left')
plt.savefig("psi_thickness.pdf",bbox_inches='tight')


""" As the plot shows the biggest difference in functions when measuring thin film thickness is around the "experimental Brewster angle" of the sample. 
The fit accuracy is therefore largely influenced by the measurements around that point."""


# test plot DELTA
angles=linspace(deg2rad(40),deg2rad(80),500)
N1plot = 1.48 + 0.0j    #SiO2
N2plot = 3.5 + 0.014j   #Si
plt.cla()
plt.plot(rad2deg(angles), rad2deg(delta(lambd, 0.3, N1plot, N2plot, angles)),'r-', label=r'd(SiO$_2$) = 0.3 nm')
plt.plot(rad2deg(angles), rad2deg(delta(lambd, 5., N1plot, N2plot, angles)),'g-', label=r'd(SiO$_2$) = 5.0 nm')
plt.plot(rad2deg(angles), rad2deg(delta(lambd, 10., N1plot, N2plot, angles)),'b-', label=r'd(SiO$_2$) = 10.0 nm')

plt.xlabel(r'$\theta_{i}$ $(^{\circ})$',size=20)
plt.ylabel(r'$\Delta$ $(^{\circ})$',size=25,rotation=90)
plt.legend(loc='lower left')
plt.savefig("delta_thickness.pdf",bbox_inches='tight')


""" FITTING THE DATA """

""" LSAT sample """

#data import 
data_import=[] 

file_name = 'LSAT.txt'
with open(file_name) as csvfile:
    readCSV = csv.reader(csvfile, delimiter='\t')
    for row in readCSV:
        data_import.append(row)

print(data_import[0]) # print header
data_import.pop(0) # .pop removes the header

# We need data formatted in two arrays: th_data=[angles of incidence] & psi_data=[psi values] 
theta_data, delta_data, psi_data, psierr_data, deltaerr_data = [], [], [], [], []

for line in data_import:
    theta_data.append(deg2rad(float(line[1]))) # be careful with deg and rad conversion!
    psi_data.append(deg2rad(float(line[2])))
    psierr_data.append(1/(deg2rad(float(line[3])))**2)
    delta_data.append(deg2rad(float(line[4])))
    deltaerr_data.append(deg2rad(float(line[5])))

psierr_data = asarray(psierr_data)
psierr_data = psierr_data/psierr_data.max()

deltaerr_data = asarray(deltaerr_data)
deltaerr_data = deltaerr_data/deltaerr_data.max()
# asarray() converts to numpy array - needed for weighting, as matehmatical operations don't work on lists
print(psierr_data)
print(deltaerr_data)

#fit SUBSTRATE -> n1=n2, k1=k2

def fcn2min(params, theta_data, data): 
    n1  = params['n1']
    k1  = params['k1']
    d   = params['d']
    
    # model - data
    residual_psi = psi(lamd, d, n1 + 1j*k1, n1 + 1j*k1, theta_data) - psi_data
    weighted = (asarray(residual_psi))**2 * psierr_data #psierr_data is already squared and inversed!
    return weighted.tolist()

# create a set of Parameters
params = Parameters()
params.add('n1',   value= 1.8, min = 1.0)
params.add('k1',   value= 0.000,  min = 0, vary=True)
params.add('d',    value= 0, vary=False)#,  min=0.02, max=0.04)


# do fit, here with LEASTSQ method
minner_leastsq = Minimizer(fcn2min, params, fcn_args=(theta_data, psi_data))
result_leastsq = minner_leastsq.minimize()
# calculate final result
final_leastsq = psi_data + result_leastsq.residual

# do fit, here with NELDER method
minner_nelder = Minimizer(fcn2min, params, fcn_args=(theta_data, psi_data))
result_nelder = minner_nelder.minimize(method= 'nelder')
# calculate final result
final_nelder = psi_data + result_nelder.residual

# write error report
print("LEASTSQ REPORT")
report_fit(result_leastsq)
print("\n \n")
print("NELDER REPORT")
report_fit(result_nelder)

# plot results
plt.cla()
plt.plot(rad2deg(theta_data), rad2deg(psi_data),'b.', label='data')
#plt.plot(rad2deg(theta_data), final_leastsq, 'r-', label='leastsq')
plt.plot(rad2deg(theta_data), rad2deg(final_nelder), 'r--', label='Nelder fit')
plt.legend(loc='lower left')

plt.xlim(44,67)

#plt.setp(fig,lw=2.5)
plt.xlabel(r'$\theta_{i}$ $(^{\circ})$',size=20)
plt.ylabel(r'$\psi$ $(^{\circ})$',size=25,rotation=90)
plt.title('LSAT substrat')
plt.savefig("psi_fit_LSAT.pdf",bbox_inches='tight')