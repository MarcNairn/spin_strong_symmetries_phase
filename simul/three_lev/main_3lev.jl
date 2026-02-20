
using DifferentialEquations, PyPlot, LinearAlgebra, FFTW, DelimitedFiles, Statistics


tr(A) = sum(diag(A))

# Gell-Mann matrices (1/2 factor)
λ1 = 1/2 * [0 1 0;
            1 0 0;
            0 0 0]

λ2 = 1/2 * [0 -im  0;
            im  0  0;
            0   0  0]

λ3 = 1/2 * [1  0 0;
            0 -1 0;
            0  0 0]

λ4 = 1/2 * [0 0 1;
            0 0 0;
            1 0 0]

λ5 = 1/2 * [0  0 -im;
            0  0   0;
            im 0   0]

λ6 = 1/2 * [0 0 0;
            0 0 1;
            0 1 0]

λ7 = 1/2 * [0  0   0;
            0  0  -im;
            0  im  0]

λ8 = (1/(2*sqrt(3))) * [1 0 0;
                        0 1 0;
                        0 0 -2]

λ_matrices = [λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8]

function su3_d_tensor(lambdas::Vector{Matrix{ComplexF64}})
    d = zeros(Float64, 8, 8, 8)
    for i in 1:8, j in 1:8, k in 1:8
        A = lambdas[i]*lambdas[j] + lambdas[j]*lambdas[i]
        d[i,j,k] = real(tr(A * lambdas[k]) / 4)
    end
    return d
end

dsu3 = su3_d_tensor(λ_matrices)

function coupled_eqs!(du, u, p, t)
    Δ1, Δ2, g, φ, Δ, κ, η, γ = p
    λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8, α_re, α_im = u

    α  = α_re + im*α_im
    αc = conj(α)
    expp = cis(φ)
    expm = conj(expp)

    d1 = -Δ1*λ2 + im*g*(λ3*(α-αc) + expm*αc/2*(-im*λ5+λ4) + expp*α/2*(-im*λ5-λ4)) + γ*(λ1*λ3 - 1/2*(λ4*λ6 + λ5*λ7))
    d2 =  Δ1*λ1 + im*g*(im*λ3*(α+αc) + expm*αc/2*(im*λ4+λ5) + expp*α/2*(im*λ4-λ5)) + γ*(λ2*λ3 + 1/2*(λ4*λ7 - λ5*λ6))
    d3 =  im*g*(-im*λ2*(α+αc) + λ1*(αc-α) + expm*αc/2*(im*λ7-λ6) + expp*α/2*(im*λ7+λ6)) + γ*(-λ1^2 - λ2^2 + 1/2*(λ6^2+λ7^2))
    d4 = (-Δ1-Δ2)*λ5 + im*g*(αc/2*(im*λ7+λ6) + α/2*(im*λ7-λ6) + expm*αc/2*(-im*λ2-λ1) + expp*α/2*(-im*λ2+λ1))
    d5 = -(-Δ1-Δ2)*λ4 + im*g*(αc/2*(-im*λ6+λ7) + α/2*(-im*λ6-λ7) + expm*αc/2*(im*λ1-λ2) + expp*α/2*(im*λ1+λ2))
    d6 = -Δ2*λ7 + im*g*(αc/2*(im*λ5-λ4) + α/2*(im*λ5+λ4) + 1/2*(expp*α - expm*αc)*(sqrt(3)*λ8 - λ3)) + γ/2*(λ1*λ4 + λ2*λ5 - λ3*λ6 + sqrt(3)*λ6*λ8)
    d7 =  Δ2*λ6 + im*g*(αc/2*(-im*λ4-λ5) + α/2*(-im*λ4+λ5) + 1/2*(expp*α + expm*αc)*(im*sqrt(3)*λ8 - im*λ3)) + γ/2*(λ1*λ5 - λ2*λ4 - λ3*λ7 + sqrt(3)*λ7*λ8)
    d8 = -im*sqrt(3)/2*g*(αc*expm*(im*λ7-λ6) + α*expp*(im*λ7+λ6)) - γ*λ8 + γ/2 * (-sqrt(3)*λ6^2 - sqrt(3)*λ7^2)

    du[1] = real(d1)
    du[2] = real(d2)
    du[3] = real(d3)
    du[4] = real(d4)
    du[5] = real(d5)
    du[6] = real(d6)
    du[7] = real(d7)
    du[8] = real(d8)

    dα = (Δ*im - κ/2)*α - im*g*((λ1 - im*λ2) + expm*(λ6 - im*λ7)) + η
    du[9]  = real(dα)
    du[10] = imag(dα)
end


g, κ = 0.1, 1.0
γ = 0.0
Δ1, Δ2 = 0.0, 0.0
Δ = 0.0
ϕ = 0.0
η = 0.5 * g

tend = 2000.0 / g
t_list = range(0.0, tend, length=10001)

# initial pure state psi -> rho0
A = [1.0, 0.0, 0.0]
O = [0.0, 1.0, 0.0]
B = [0.0, 0.0, 1.0]
a = 1/sqrt(2) # free to choose 
b = -1/sqrt(2) # free to choose as long as normalization is conserved (|b|≤|a|)
c = sqrt(max(0.0, 1 - a^2 - b^2)) # normalization
psi = a*A + b*B + c*O
rho0 = psi * psi'
rho0 = rho0 / tr(rho0)

u0 = zeros(Float64, 10)
for i in 1:8
    u0[i] = real(tr(rho0 * λ_matrices[i]))
end

params = (Δ1, Δ2, g, ϕ, Δ, κ, η, γ)
prob = ODEProblem(coupled_eqs!, u0, (0.0, tend), params)
sol = solve(prob, Tsit5(), saveat=collect(t_list), reltol=1e-9, abstol=1e-9)

λ = sol[1:8, :]
t = sol.t

N0 = 1/3 .- sol[3,:] .+ 1/sqrt(3).*sol[8,:]
N1 = 1/3 .+ sol[3,:] .+ 1/sqrt(3).*sol[8,:]
N2 = 1/3 .- 2/sqrt(3).*sol[8,:]
C2 = sum(λ .^ 2, dims=1)[:]

# cubic Casimir
Ntsteps = length(t)
C3 = zeros(Float64, Ntsteps)
@inbounds for ti in 1:Ntsteps
    val = 0.0
    for i in 1:8, j in 1:8, k in 1:8
        val += dsu3[i,j,k] * λ[i,ti] * λ[j,ti] * λ[k,ti]
    end
    C3[ti] = val
end

# save
#mkpath("simul")
#writedlm("simul/N0.txt", N0)
#writedlm("simul/N1.txt", N1)
#writedlm("simul/N2.txt", N2)

using NPZ

# grid phase diagram
a_vals = range(-1.0, 1.0, length=151)
b_vals = range(-1.0, 1.0, length=151)


tend = 100.0 / g
t_list = range(0.0, tend, length=1001)
steady_idx = collect(div(length(t_list), 2):length(t_list))

eta_ratios = 0.8:-0.01:0.01 #slices

data_slices = Dict{Float64, Array{Float64,2}}()

for (progress, η_ratio) in enumerate(eta_ratios)
    η_val = η_ratio * g
    params = (0.0, 0.0, g, 0.0, Δ, κ, η_val, γ)
    peak_amp = zeros(Float64, length(a_vals), length(b_vals))

    for (i, a) in enumerate(a_vals)
        for (j, b) in enumerate(b_vals)
            if a^2 + b^2 ≤ 0.995
                c = sqrt(max(0.0, 1.0 - a^2 - b^2))
                ψ = [a, c, b]
                ρ0 = ψ * ψ'
                ρ0 /= tr(ρ0)

                u0_scan = zeros(Float64, 10)
                for k in 1:8
                    u0_scan[k] = real(tr(ρ0 * λ_matrices[k]))
                end

                prob = ODEProblem(coupled_eqs!, u0_scan, (0.0, tend), params)
                sol = solve(prob, Tsit5(); saveat=collect(t_list), reltol=1e-6, abstol=1e-6)

                N2_ss = @. 1/3 - (2/√3) * sol[8, steady_idx]
                N2_ss .-= mean(N2_ss)
                fft_N2 = fft(N2_ss)
                if length(fft_N2) > 1
                    amp = maximum(abs.(fft_N2[2:end])) #order parameter
                else
                    amp = 0.0
                end
                peak_amp[i, j] = amp
            else
                peak_amp[i, j] = 0.0
            end
        end
    end

    maxval = maximum(peak_amp)
    if maxval > 0.0
        normalized_amp = peak_amp ./ maxval
    else
        normalized_amp = peak_amp
    end

    data_slices[η_ratio] = normalized_amp
    println("Completed $(progress)/$(length(eta_ratios)) (η/g = $(round(η_ratio, digits=3)))")
end

save_dict = Dict{String, Any}()
save_dict["a_vals"] = collect(a_vals)
save_dict["b_vals"] = collect(b_vals)
save_dict["eta_ratios"] = collect(eta_ratios)
for η_ratio in eta_ratios
    key = "data_slice_$(η_ratio)"
    save_dict[key] = data_slices[η_ratio]
end

npzwrite("ab_eta_slices.npz", save_dict)
```
