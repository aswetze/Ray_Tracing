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

μ2, dμ2_dr = index_refraction_no_b(6677000, 5e6, ["QP", 5.e11, 300, 50])
println(μ2, " ", dμ2_dr)
