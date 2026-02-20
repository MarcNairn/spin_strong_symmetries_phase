using DifferentialEquations
using LaTeXStrings


default(
    fontfamily = "Computer Modern", # TeX's default serif font
    tickfont = font(12, "Computer Modern"),  # Tick labels font
    guidefont = font(14, "Computer Modern"), # Axis labels font
    legendfont = font(12, "Computer Modern"), # Legend font
    titlefont = font(16, "Computer Modern")  # Title font
)

# Define parameters 
δ = 0.5*Δ
Δ = 1.0
κ = 1.0
g = 3.0              
chi_R = 1/(im*Δ - κ/2)
θ = atan(2Δ/κ)             
ϕ = π/3               
Δ0 = 1.0 + δ       
Δ1 = 1.0 - δ               
Γ_T = 1.0              
Γ_up = 0.60*Γ_T          
Γ_down = Γ_T - Γ_up           

# Compute interaction parameters
χ = g^2 * abs(chi_R) * exp(im * θ)
ξ₊ = g^2 * abs(chi_R) * exp(im * θ) * exp(im * ϕ)
ξ₋ = g^2 * abs(chi_R) * exp(im * θ) * exp(-im * ϕ)

# Define the system of ODEs
function system!(du, u, p, t)
    s0m, s1m, s0z, s1z = u  # State variables
    
    # time evo
    du[1] = - (im * Δ0 + Γ_T / 2) * s0m + χ * s0m * s0z + ξ₊ * s1m * s0z
    du[2] = - (im * Δ1 + Γ_T / 2) * s1m + ξ₋ * s0m * s1z + χ * s1m * s1z
    du[3] = -4 * real(conj(χ) * conj(s0m) * s0m + conj(ξ₊) * conj(s1m) * s0m) - Γ_T * s0z + Γ_up - Γ_down
    du[4] = -4 * real(conj(ξ₋) * conj(s0m) * s1m + conj(χ) * conj(s1m) * s1m) - Γ_T * s1z + Γ_up - Γ_down
end

# Initial conditions (modify as needed)
s0m₀ = 0 
s1m₀ = 0
s0z₀ = 1 - Γ_up/Γ_T
s1z₀ = 1 - Γ_up/Γ_T
u₀ = [s0m₀, s1m₀, s0z₀, s1z₀]

# Time span
tend= 10/κ
tspan = (0.0, tend)

# Solve the system
prob = ODEProblem(system!, u₀, tspan; saveat=tend/100)
sol = solve(prob)

sz0 = sol[3,:]
sz1 = sol[4, :]

sm0 = sol[1,:]
sm1 = sol[2, :]

plot(sol.t, [sm0,sm1,sz0,sz1], linewidth=5, xlabel = L"t\kappa", ylabel = L"\sigma^\alpha", label= [L"\sigma^-_0" L"\sigma^-_1" L"\sigma^z_0" L"\sigma^z_1" ])

#sol_ss = SteadyStateProblem(prob)

    