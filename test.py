import numpy
import matplotlib.pylab as plt
import scipy
from scipy import integrate


def QPNe(r, Nm, Hm, Ym):
    # Nm in m^{-3}, Hm in km, Ym in km
    # croft and hoogasian ionosphere - good for checking
    # https://doi.org/10.1002/rds19683169
    HmMeters = Hm * 1000.
    YmMeters = Ym * 1000.
    Rearth = 6371. * 1000.  # m
    rMax = (HmMeters + Rearth)
    rBase = rMax - YmMeters
    rTop = rMax + YmMeters

    if (r > rBase) and (r < rTop):
        a = (r - rMax) / (YmMeters)
        Ne = Nm * (1. - (a * a))
        dNedr = -2. * Nm * a / YmMeters
    else:
        Ne = 1.e-31
        dNedr = 1.e-31
    return Ne, dNedr


# set some of these things in a class or at initialization
def IndexRefractionNoB(r, frequency, IonoParams):
    # frequency in Hz
    # mu squared
    #     print(IonoParams)
    if IonoParams[0] == 'QP':
        Nm = float(IonoParams[1])
        Hm = float(IonoParams[2])
        Ym = float(IonoParams[3])
        Ne, dNedr = QPNe(r, Nm, Hm, Ym)
    fePlasma = 8.98e3 * numpy.sqrt(Ne / 1.e6)
    X = fePlasma / frequency  # to be consistent with appleton hartree equation
    mu2 = 1. - X * X
    dXdr = X * dNedr / (2. * Ne)
    dmu2dr = -2. * X * dXdr

    return mu2, dmu2dr

Nm = 5.e11
Ym = 50.
Hm = 300.
frequency = 5.e6

rTest = numpy.arange(6371.,6700.,1.)*1000.

for ii in rTest:
    Ne,dNedr = QPNe(ii, Nm,Hm,Ym)
    print(ii, Ne, dNedr)

IonoParams = numpy.array(['QP',Nm,Hm,Ym])

for ir in rTest:
    mu2,dmu2dr = IndexRefractionNoB(ir,frequency, IonoParams)
    print(ir, mu2,dmu2dr)


def Lagrangian(t, x, args):
    #     print(x)
    r = x[0]
    theta = x[1]
    Q = x[2]
    frequency = args[0]
    IonoParams = args[1]
    y = numpy.zeros(3)
    # call index of refraction
    mu2, dmu2dr = IndexRefractionNoB(r, frequency, IonoParams)

    # P1 is the group delay
    dQdP1 = dmu2dr / 2. + (mu2 - Q * Q) / r  # equation 10
    drdP1 = Q
    dthetadP1 = numpy.sqrt(mu2 - Q * Q) / r
    #     print(y)
    y[0] = drdP1
    y[1] = dthetadP1
    y[2] = dQdP1
    return y

frequency = 8.e6
Nm = 5.e11
Ym = 50.
Hm = 300.
r = 6371.e3
IonoParams = numpy.array(['QP',Nm,Hm,Ym])
IndexRefractionNoB(r,frequency, IonoParams)

# set initial conditions
Elevation = 45.
Qinitial = numpy.sin(numpy.deg2rad(Elevation))
rInitial = 6371.e3
thetaInitial = 0.
x0 = numpy.array([rInitial,thetaInitial,Qinitial])

print

# set parameters
frequency = 8.e6
Nm = 5.e11
Ym = 50.
Hm = 300.
IonoParams = numpy.array(['QP',Nm,Hm,Ym])

print(Lagrangian(0.,x0, (frequency, IonoParams)))
# set initial vector

t0 = 0.
# initalize ode integration with adapative step size
# little help from other program...
# need to help each time
rdopri = scipy.integrate.ode(Lagrangian).set_integrator('dopri5',\
                    atol=1e-6, first_step=10e3,max_step=10e3,dfactor=0.1,nsteps=1,)
rdopri.set_f_params((frequency,IonoParams))
rdopri.set_initial_value(x0,t0)

x = numpy.zeros([3,1000])*numpy.nan
t1 = numpy.zeros(1000)*numpy.nan
k = 0
tIntLimit = 8000.*1e3
t = t0
k=0
while t < tIntLimit:
    tx = rdopri.integrate(tIntLimit,step=True)
    print(tx)
    t = rdopri.t
    print(t)
    t1[k] = rdopri.t
    x[:,k] = tx
    k+=1
    if tx[0] < 6371e3:
        break
    if tx[0] > 6871e3:
        break
    if k > 998:
        break

dist = x[1,:]*6371e3
height = x[0,:] - 6371e3
plt.figure(dpi=300)
plt.plot(dist,height, 'k.')
plt.savefig("plots_test.pdf")
