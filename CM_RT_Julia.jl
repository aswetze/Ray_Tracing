using DifferentialEquations

#Checked for correctness against Coleman_RT_Python.ipynb
function qpne(r, nm, hm, ym)
    hm_meters = hm * 10^3
    ym_meters = ym * 10^3
    r_earth = 6371 * 10^3
    r_max = hm_meters + r_earth
    r_base = r_max - ym_meters
    r_top = r_max + ym_meters

    if (r > r_base && r < r_top)
        a = (r - r_max) / ym_meters
        ne = nm * (1. - a^2)
        dne_dr = -2. * nm * a / ym_meters
    else
        ne = 1. * 10^-31
        dne_dr = 1. * 10^-31
    end
    return ne, dne_dr
end

#Checked for correctness against Coleman_RT_Python.ipynb
function index_refraction_no_b(r, freq, iono_params)
    nm = iono_params[2]
    hm = iono_params[3]
    ym = iono_params[4]
    ne, dne_dr = qpne(r, nm, hm, ym)


    fe_plasma = 8.98e3 * sqrt(ne / 1.e6)
    x = fe_plasma / freq
    mu2 = 1. - x * x
    dx_dr = x * dne_dr / (2. * ne)
    dmu2_dr = -2. * x * dx_dr
    return mu2, dmu2_dr
end

#Test arrays
rTest = collect(6371000.:1000:6699000)
IonoParams = ["QP", 5.e11, 300., 50.]

#Tests function qpne
for i in rTest
    Ne, dNedr = qpne(i, 5.e11, 300., 50.)
    if Ne == 1.0000000000000034e-31
        Ne = "Default"
    end
    if dNedr == 1.0000000000000034e-31
        dNedr = "Default"
    end
    #println(i, " ", Ne, " ", dNedr)
end

#Tests function index_refraction_no_b
for i in rTest
    mu2, dmu2dr = index_refraction_no_b(i, 5.e6, IonoParams)
    if mu2 == 1.0
        mu2 = "Default"
    end
    #println(i, " ", mu2, " ", dmu2dr)
end

function raytrace!(du, u, p, t)
    q = u[1]
    r = u[2]
    θ = u[3]

    #μ^2 = -0.581709 + (-1.68825*10^-14 + 1.26537*10^-21 r) r
    #dμ2/dr = -0.00844126 + 1.26537*10^-9 r
    if r > 6621000 && r < 6721000
        du[1] = (-.00844126+1.26537*10^-9*r)/2 + ((-.0581709 + (-1.68825*10^-14+1.26537*10^-21*r)*r)-q^2)/r
        du[3] = sqrt(abs((-.0581709 + (-1.68825*10^-14+1.26537*10^-21*r)*r)-q^2))/r
    else
        du[1] = (-3.225616*10^-43)/2 + (1-q^2)/r
        du[3] = sqrt(abs(1 - q^2)) / r
    end
    du[2] = q
end

#u0 = Float64[q0, r0, Θ0]
u0 = Float64[sin(π/4); 6371.e3; 0.0]
tspan = (0.0f0, 748000.0f0)
prob = ODEProblem(raytrace!, u0, tspan, dtmax=100)
sol = solve(prob)
#plot(sol)
# println(sol[1,:])
# println("sol2 = ")
# println(sol[2,:])
# println("sol3 = ")
# println(sol[3,:])
# println("sol.t")
# println(sol.t)
# println(sol[2,:].-6371e3)
plot(sol.t .* 6371, sol[2,:].-6371.e3)
