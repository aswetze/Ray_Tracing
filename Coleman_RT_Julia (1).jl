using OrdinaryDiffEq, DiffEqGPU, Plots, ODEInterfaceDiffEq, ODEInterfaceDiffEq

function qpne(r, nm, hm, ym)
    hm_meters = hm * 10^3
    ym_meters = ym * 10^3
    r_earth = 6371 * 10^3
    r_max = hm_meters + r_earth
    r_base = r_max - ym_meters
    r_top = r_max + ym_meters

    if (r > r_base && r < r_top)
        a = (r - r_max) / ym_meters
        ne = nm * (1 - a^2)
        dne_dr = -2 * nm * a / ym_meters
    else
        ne = 1 * 10^-31
        dne_dr = ne
    end
    return ne, dne_dr
end

function index_refraction_no_b(r, freq, iono_params)
    nm = iono_params[2]
    hm = iono_params[3]
    ym = iono_params[4]
    ne, dne_dr = qpne(r, nm, hm, ym)

    fe_plasma = 8.98e3 * sqrt(ne / 1.e6)
    x = fe_plasma / freq
    μ2 = 1. - x * x
    dx_dr = x * dne_dr / (2. * ne)
    dμ2_dr = -2. * x * dx_dr
    return μ2, dμ2_dr
end

function raytrace!(du, u, p, t)
    q = u[1]
    r = u[2]
    θ = u[3]
    μ2, dμ2_dr = index_refraction_no_b(r, 8.e6, ["QP", 5.e10, 300, 50])

    du[1] = dμ2_dr/2 + (μ2-q^2)/r
    du[2] = q
    du[3] = sqrt(μ2-q^2)/r
end


# Main
#******************************************************************************
r0 = 6371.e3
Q0 = sin(π/4)
θ0 = 0.0

u0 = Float64[Q0; r0; θ0]
tspan = (0.0f0,12000.f0)
prob = ODEProblem(raytrace!, u0, tspan)
@time sol = solve(prob, Tsit5())
#println(sol[1,:])
#println("sol2 = ")
println(sol[2,:] .- 6371e3)
#println("sol3 = ")
#println(sol[3,:])
println("sol.t")
println(sol.t)
plot(sol.t .* 6371e3, sol[2,:] .- 6371e3)
