using Plots, DifferentialEquations
function raytrace!(du, u, p, t)
    q = u[1]
    r = u[2]
    θ = u[3]

    μa = 4.94284*10^-22
    μb = -6.59474*10^-15
    μc = 0.0114319
    dμa = 7.90854*10^-10
    dμb = -0.00527579

    if r > 6621000 && r < 6721000
        du[1] = (dμa*r + dμb) /2 + ((μa*r^2 + μb*r + μc)-q^2)/r
        du[3] = sqrt(abs((μa*r^2 + μb*r + μc)-q^2))/r
    else
        du[1] = (-3.225616*10^-43)/2 + (1-q^2)/r
        du[3] = sqrt(abs(1 - q^2)) / r
    end
    du[2] = q
end

u0 = Float64[sin(π/4); 6371.e3; 0.0]
tspan = (0.0f0, 790000.0f0)
prob = ODEProblem(raytrace!, u0, tspan, dtmax=1000)
sol = solve(prob)
plot(sol.t *6371, sol[2,:].-6371.e3)
