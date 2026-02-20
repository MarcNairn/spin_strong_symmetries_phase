using QuantumCumulants
using ModelingToolkit, OrdinaryDiffEq
using Plots
using LaTeXStrings
using PGFPlotsX
default(
    fontfamily = "Palatino", 
    tickfont = font(12, "Palatino"),  # Tick labels font
    guidefont = font(14, "Palatino"), # Axis labels font
    legendfont = font(12),#, "Computer Modern"), # Legend font
    titlefont = font(16)#, "Computer Modern")  # Title font
)

# Define parameters
κ, g, Δ, Δ1, Δ2, ϕ, η, Γ= cnumbers("κ g Δ Δ1 Δ2 ϕ η Γ")

N=2

# Define hilbert space
hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom,i),2) for i=1:N]...)
h = hf ⊗ ha

# Define the fundamental operators
a = Destroy(h,:a)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k+1)

# Define the Hamiltonian
function pickH(name::String)
    if name=="TC"
        H = -Δ * a' * a - Δ1 * σ(2,2,1) - Δ2*σ(2,2,2) + sum(g*(exp(-im*ϕ*(i-1))*a'*σ(1,2,i) + exp(im*ϕ*(i-1))*a*σ(2,1,i)) for i=1:N) + im*η*(a'-a)
    elseif name=="Dicke"
        H = -Δ * a' * a + (Δ1 - g*(a'+a))*σ(2,2,1) + (Δ2- g*(exp(-im*ϕ)*a'+exp(im*ϕ)*a))*σ(2,2,2) + sum(g*(exp(-im*ϕ*(i-1))*a'*σ(1,2,i) + exp(im*ϕ*(i-1))*a*σ(2,1,i) +exp(-im*ϕ*(i-1))*a'*σ(2,1,i) + exp(im*ϕ*(i-1))*a*σ(1,2,i)) for i=1:N)
    end
return H
end


H = pickH("TC")
# Collapse operators
J = [a]#, σ(2,1,1), σ(2,1,2), σ(1,2,1), σ(1,2,2)]
rates = [κ]#, 0.75*Γ, 0.75*Γ, 0.25/2*Γ, 0.25/2*Γ]

# Define the operators of interest

sx1 = σ(1,2,1) + σ(2,1,1)
sx2 = σ(1,2,2) + σ(2,1,2)
sy1 = im*(σ(1,2,1) - σ(2,1,1))
sy2 = im*(σ(1,2,2) - σ(2,1,2))
sz1 = σ(1,1,1) - σ(2,2,1)
sz2 = σ(1,1,2) - σ(2,2,2)


# ops = [a,a',sx1, sx2, sy1, sy2, sz1, sz2]
ops = [a, a',  σ(1,2,1), σ(1,2,2), σ(2,2,1), σ(2,2,2)]

# Derive equation for average photon number
eqs = meanfield(ops,H,J;rates=rates,order=1)
#eqs_expanded = cumulant_expansion(eqs,2)

me_comp = complete(eqs) #automatically complete the system

# Build an ODESystem out of the MeanfieldEquations
@named sys = ODESystem(me_comp)


# initial state'
u0 = zeros(ComplexF64, length(me_comp))
u0[end-1] = -1
u0[end] = -1


δ = 0;
ω0 = 1.0

κn= 1*ω0
gn = 0.1*ω0
Δn = 0;
Δ1n = 0*ω0 - δ
Δ2n = 0*ω0 + δ
ηn = 0.2*ω0
ϕn = 0



# list of parameters
ps = (Δ, g, κ, ϕ, Δ1, Δ2, η)
p0 = ps .=> (Δn, gn, κn, ϕn, Δ1n, Δ2n, ηn)
tend = 3000/κn

prob = ODEProblem(sys,u0,(0.0,tend),p0)# Solve the initial ODE problem
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)


""" PLOTS """

s1z = [real.(2 * sol[σ(2, 2, 1)] .- 1)]
s2z = [real.(2 * sol[σ(2, 2, 2)] .- 1)]
s1x = [real.(sol[σ(1, 2, 1)] + sol[σ(2, 1, 1)])]
s2x = [real.(sol[σ(1, 2, 2)] + sol[σ(2, 1, 2)])]
s1y = [imag.(sol[σ(1, 2, 1)] - sol[σ(2, 1, 1)])]
s2y = [imag.(sol[σ(1, 2, 2)] - sol[σ(2, 1, 2)])]
a_phase = [angle.(sol[a])]
a_dag_a = [real.(sol[a'] .* sol[a])]


plot(sol.t, [s1z, s2z], xlabel = L"t\kappa", ylabel = L"\sigma^α", linewidth = 5, legend=:topright, color = [:teal :steelblue], alpha=0.7, linestyle= [:solid :solid])

# Extract the final state
final_state = sol.u[end]





# Update parameters with new ϕ value
p0_new = ps .=> (Δn, gn, κn,-1*ϕn, Δ1n, Δ2n, ηn)

# Solve the new ODE problem with the final state as the initial state
prob = ODEProblem(sys,final_state,(0.0,tend),p0_new)
sol_new = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Concatenate time and solution arrays
t_combined = [sol.t; sol_new.t .+ tend]
a_phase = [angle.(sol[a]); angle.(sol_new[a])]
a_dag_a = [real.(sol[a'] .* sol[a]); real.(sol_new[a'] .* sol_new[a])]
sigmaz1_combined = [real.(2 * sol[σ(2, 2, 1)] .- 1); real.(2 * sol_new[σ(2, 2, 1)] .- 1)]
sigmaz2_combined = [real.(2 * sol[σ(2, 2, 2)] .- 1); real.(2 * sol_new[σ(2, 2, 2)] .- 1)]
sigmax1_combined = [real.(sol[σ(1, 2, 1)] + sol[σ(2, 1, 1)]); real.(sol_new[σ(1, 2, 1)] + sol_new[σ(2, 1, 1)])]
sigmax2_combined = [real.(sol[σ(1, 2, 2)] + sol[σ(2, 1, 2)]); real.(sol_new[σ(1, 2, 2)] + sol_new[σ(2, 1, 2)])]
sigmay1_combined = [imag.(sol[σ(1, 2, 1)] - sol[σ(2, 1, 1)]); imag.(sol_new[σ(1, 2, 1)] - sol_new[σ(2, 1, 1)])]
sigmay2_combined = [imag.(sol[σ(1, 2, 2)] - sol[σ(2, 1, 2)]); imag.(sol_new[σ(1, 2, 2)] - sol_new[σ(2, 1, 2)])]
# Plot
plot(t_combined, [sigmaz1_combined, sigmaz2_combined, sigmax1_combined, sigmax2_combined, sigmay1_combined, sigmay2_combined], xlabel = L"t\kappa", ylabel = L"\sigma^α", label = [L"\sigma^z_1" L"\sigma^z_2" L"\sigma^x_1" L"\sigma^x_2" L"\sigma^y_1" L"\sigma^y_2"], linewidth = 5, legend=:topright, color = [:seagreen :sienna :darkcyan :darkblue :black :maroon], alpha=0.7, linestyle= [:solid :solid :dot :dot :dot :dot])
plot(t_combined, [sigmaz1_combined, sigmaz2_combined], xlabel = L"t\kappa", ylabel = L"\sigma^α", label = [L"\sigma^z_1" L"\sigma^z_2"], linewidth = 5, legend=:topright, color = [:seagreen :sienna], alpha=0.7, linestyle= [:solid :solid])

plot(t_combined, a_phase/π, xlabel = L"t\kappa", ylabel = L"\phi_l/\pi", label = "Cavity phase", linewidth = 5, legend=:topright, color = :black, alpha=0.7, linestyle= :solid)
plot(t_combined, a_dag_a, xlabel = L"t\kappa", ylabel = L"\langle a^\dagger a \rangle", label = L"\langle a^\dagger a \rangle", linewidth = 5, legend=:topright, color = :black, alpha=0.7, linestyle= :solid)