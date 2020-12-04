import numpy
import matplotlib.pylab as plt
import scipy.integrate
import sys


def qpne(r, nm, hm, ym):
    # Nm in m^{-3}, Hm in km, Ym in km
    # croft and hoogasian ionosphere - good for checking
    # https://doi.org/10.1002/rds19683169
    hm_meters = hm * 1000.
    ym_meters = ym * 1000.
    r_earth = 6371. * 1000.  # m
    r_max = (hm_meters + r_earth)
    r_base = r_max - ym_meters
    r_top = r_max + ym_meters

    if (r > r_base) and (r < r_top):
        a = (r - r_max) / ym_meters
        ne = nm * (1. - (a * a))
        dne_dr = -2. * nm * a / ym_meters
    else:
        ne = 1.e-31
        dne_dr = 1.e-31
    return ne, dne_dr


def index_refraction_no_b(r, freq, iono_params):
    # frequency in Hz
    # mu squared
    # print(IonoParams)
    if iono_params[0] == 'QP':
        nm = float(iono_params[1])
        hm = float(iono_params[2])
        ym = float(iono_params[3])
        ne, dne_dr = qpne(r, nm, hm, ym)

        fe_plasma = 8.98e3 * numpy.sqrt(ne / 1.e6)
        x = fe_plasma / freq  # to be consistent with appleton hartree equation
        mu2 = 1. - x * x
        dx_dr = x * dne_dr / (2. * ne)
        dmu2_dr = -2. * x * dx_dr
        return mu2, dmu2_dr
    else:
        print("************************** Index Refraction Failed **************************")
        return "nan", "nan"


def lagrangian(t, x, args):
    #     print(x)
    r = x[0]
    theta = x[1]
    q = x[2]
    freq = args[0]
    iono_params = args[1]
    y = numpy.zeros(3)
    # call index of refraction
    mu2, dmu2_dr = index_refraction_no_b(r, freq, iono_params)

    # P1 is the group delay
    dq_dp1 = dmu2_dr / 2. + (mu2 - q * q) / r  # equation 10
    dr_dp1 = q

    dtheta_dp1 = numpy.sqrt(mu2 - q * q) / r
    #     print(y)
    y[0] = dr_dp1
    y[1] = dtheta_dp1
    y[2] = dq_dp1
    return y


# ********************* MAIN *********************

# set initial conditions
Elevation = 45.
Qinitial = numpy.sin(numpy.deg2rad(Elevation))
rInitial = 6371.e3
thetaInitial = 0.
x0 = numpy.array([rInitial, thetaInitial, Qinitial])

# set parameters
frequency = 8.e6
Nm = 5.e11
Ym = 50.
Hm = 300.
IonoParams = numpy.array(['QP', Nm, Hm, Ym])

t0 = 0.
# initalize ode integration with adapative step size
# little help from other program...
# need to help each time
rdopri = scipy.integrate.ode(lagrangian).set_integrator('dopri5', atol=1e-6, first_step=10e3, max_step=10e3, dfactor=0.1, nsteps=1, )
rdopri.set_f_params((frequency, IonoParams))
rdopri.set_initial_value(x0, t0)

X = numpy.zeros([3, 1000])*numpy.nan
t1 = numpy.zeros(1000)*numpy.nan
K = 0
tIntLimit = 8000.*1e3
T = t0


while T < tIntLimit:
    tx = rdopri.integrate(tIntLimit, step=True)
    T = rdopri.t
    t1[K] = rdopri.t
    X[:, K] = tx
    K += 1
    if tx[0] < 6371e3:
        break
    if tx[0] > 6871e3:
        break
    if K > 998:
        break


# ********************* OUTPUT *********************

dist = X[1, :]*6371e3
height = X[0, :] - 6371e3
plt.figure(dpi=300)
plt.plot(dist, height, 'k.')
plt.savefig("plots.pdf")

# ********************* END *********************
sys.stdout.close()
