using DifferentialEquations
using PyPlot
using LaTeXStrings
using LinearAlgebra
using Statistics
using ProgressMeter
using PyCall
using FFTW
using DelimitedFiles

PyPlot.pygui(false)
println(PyCall.python)

# ENV["PYTHON"] = "C:\\Users\\admin\\Anaconda3\\envs\\data_vis\\python.exe"

LogNorm = pyimport("matplotlib.colors").LogNorm
cmap = pyimport("cmap")

lajolla = cmap.Colormap("crameri:lajolla").to_mpl();
balance = cmap.Colormap("cmocean:balance").to_mpl();
amethyst = cmap.Colormap("cmasher:amethyst").to_mpl();
cet_cbl = cmap.Colormap("colorcet:CET_CBL3").to_mpl();

# Configure LaTeX rendering and fonts
PyPlot.matplotlib[:rc]("text", usetex=true)
PyPlot.matplotlib[:rc]("font", family="serif", serif=["mathpazo"], size=18)  # Base font size
PyPlot.matplotlib[:rc]("axes", titlesize=22)             # Axis title
PyPlot.matplotlib[:rc]("axes", labelsize=20)             # Axis labels
PyPlot.matplotlib[:rc]("xtick", labelsize=18)            # X-ticks
PyPlot.matplotlib[:rc]("ytick", labelsize=18)            # Y-ticks
PyPlot.matplotlib[:rc]("legend", fontsize=18)            # Legend
PyPlot.matplotlib[:rc]("figure", titlesize=24)           # Figure title
PyPlot.svg(true)
# LaTeX preamble packages
PyPlot.matplotlib[:rc]("text.latex", preamble="\\usepackage{amsmath}\\usepackage{amsfonts}\\usepackage{amssymb}\\usepackage{lmodern}")


""" 
#############################################################

                    USEFUL FUNCTIONS BELOW

#############################################################
"""


function compute_time_average(sol::ODESolution{}, variable_index::Int)
    """  
    ###################
    Compute time average of an ODE sol object using portion of data, requires indexed variables 
    ################### 
    """
    t = sol.t
    t_mid = 0.6*t[end]
    indices = findall(x -> x >= t_mid, t)
    vals = sol[variable_index, indices]
    return mean(vals), var(vals)
end


function compute_time_average(t::Vector{Float64}, series::Vector{Float64})
    """
    ################### 
    Compute time average of a generic series using portion of data, no indexing required
    ################### 
    """
    t_mid = 0.6 * maximum(t)
    indices = findall(x -> x >= t_mid, t)
    return mean(series[indices]), var(series[indices])
end

"""
################################
DYNAMICS 
################################
"""

function spin_cav_ode!(du, u, p, t)
    s1x, s2x, s1y, s2y, s1z, s2z, α_r, α_i = u
    Δ1, Δ2, Δ, g, φ, κ, η = p

    # Equations for spin 1
    du[1] = Δ1 * s1y - 2*g * s1z * α_i            # ds1x/dt
    du[3] = -Δ1 * s1x - 2*g * s1z * α_r            # ds1y/dt
    du[5] = 2g * (s1x * α_i + s1y * α_r)          # ds1z/dt

    # Equations for spin 2
    du[2] = Δ2 * s2y - 2g * s2z * (sin(φ)*α_r + cos(φ)*α_i)  # ds2x/dt
    du[4] = -Δ2 * s2x   - 2g * s2z * (cos(φ)*α_r - sin(φ)*α_i) # ds2y/dt
    du[6] = 2g * (s2x*(cos(φ)*α_i + sin(φ)*α_r) + s2y*(cos(φ)*α_r - sin(φ)*α_i)) # ds2z/dt

    # Cavity field equations
    du[7] = -Δ*α_i - (κ/2)*α_r - g/2*(s1y +  cos(φ)*s2y + sin(φ)*s2x) + η  # dα_r/dt
    du[8] = Δ*α_r - (κ/2)*α_i - g/2*(s1x - sin(φ)*s2y + cos(φ)*s2x) # dα_i/dt
    return nothing
end


function single_spin_cav_ode!(du, u, p, t)
    s1x, s1y, s1z, α_r, α_i = u
    Δ1, Δ2, Δ, g, φ, κ, η = p

    g = 2*g

    # Equations for spin 1
    du[1] = Δ1 * s1y - 2*g * s1z * α_i            # ds1x/dt
    du[2] = -Δ1 * s1x - 2*g * s1z * α_r            # ds1y/dt
    du[3] = 2g * (s1x * α_i + s1y * α_r)          # ds1z/dt

    # Cavity field equations
    du[4] = -Δ*α_i - (κ/2)*α_r - g/2*s1y + η  # dα_r/dt
    du[5] = Δ*α_r - (κ/2)*α_i - g/2*s1x  # dα_i/dt
    return nothing
end

function spinonly_ode!(du, u, p, t)
    s1x, s2x, s1y, s2y, s1z, s2z = u
    Δ1, Δ2, Δ, g, φ, κ, η = p

    #Check normalization of "g"


    # Derived parameters
    K = κ / 2
    Δ_plus = Δ * cos(φ) + K * sin(φ)
    Δ_minus = Δ * cos(φ) - K * sin(φ)
    K_plus = K * cos(φ) + Δ * sin(φ)
    K_minus = K * cos(φ) - Δ * sin(φ)
    𝒥 = g^2 / (Δ^2 + K^2)  # 𝒥 = g²/(Δ² + (κ/2)²)

    # Equations for s1x, s2x, s1y, s2y
    du[1] = Δ1 * s1y + s1z * 𝒥 * (Δ * s1y + K * s1x + K_plus * s2x + Δ_plus * s2y - 2 * K * η / g)
    du[2] = Δ2 * s2y + s2z * 𝒥 * (Δ * s2y + K * s2x + K_minus * s1x + Δ_plus * s1y - 2 * K_minus * η / g)
    du[3] = -Δ1 * s1x + s1z * 𝒥 * (-Δ * s1x + K * s1y - Δ_minus * s2x + K_plus * s2y + 2 * Δ * η / g)
    du[4] = -Δ2 * s2x + s2z * 𝒥 * (-Δ * s2x + K * s2y - Δ_plus * s1x + K_minus * s1y + 2 * Δ_plus * η / g)

    # Equations for s1z, s2z
    du[5] = -𝒥 * (
        -2 * η / g * (K * s1x - Δ * s1y) +
        K * (s1x^2 + s1y^2) +
        K_plus * s1x * s2x +
        Δ_plus * s1x * s2y -
        Δ_minus * s1y * s2x +
        K_plus * s1y * s2y
    )

    du[6] = -𝒥 * (
        -2 * η / g * (K_minus * s2x - Δ_plus * s2y) +
        K * (s2x^2 + s2y^2) +
        K_minus * s2x * s1x +
        Δ_plus * s2x * s1y -
        Δ_plus * s2y * s1x +
        K_minus * s2y * s1y
    )

    return nothing
end


function single_spinonly_ode!(du, u, p, t)
    s1x, s1y,s1z = u
    Δ1, Δ2, Δ, g, φ, κ, η = p #Keep with extra parameters just for simplicity (wont be used)

    # Derived parameters
    K = κ / 2

    𝒥 = g^2 / (Δ^2 + K^2)  # 𝒥 = g²/(Δ² + (κ/2)²)

    # Equations for s1x, s2x, s1y, s2y
    du[1] = Δ1 * s1y + s1z * 𝒥 * (Δ * s1y + K * s1x - 2 * K * η / g)
    du[2] = -Δ1 * s1x + s1z * 𝒥 * (-Δ * s1x + K * s1y + 2 * Δ * η / g)
    # Equations for s1z, s2z
    du[3] = -𝒥 * (-2 * η / g * (K * s1x - Δ * s1y) +
        K * (s1x^2 + s1y^2)
    )
    return nothing
end



""" 
##############################
SPIN OBSERVABLES  
##############################
"""

function plot_sz_phase_diagram(φ_vals, η_vals, s1z_avg, s2z_avg, s1z_vrc, s2z_vrc)
    # Create figure with subplots
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8.5), dpi=300,
                            constrained_layout=true, sharex="col")

    # Calculate extent for imshow [xmin, xmax, ymin, ymax]
    extent = [φ_vals[1], φ_vals[end], η_vals[1]/g, η_vals[end]/g]

    # Top row plots
    ax1_top, ax2_top = axs[1, 1], axs[1, 2]

    # Plot s1z_avg
    im1 = ax1_top.imshow(s1z_avg, cmap=balance, vmin=-1, vmax=1,
                        extent=extent, aspect="auto", origin="lower")
    ax1_top.set_ylabel(L"\eta/g")
    # ax1_top.set_xlabel(L"\varphi")
    ax1_top.set_xticks([-pi/2, 0, pi/2])
    ax1_top.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax1_top.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax1_top.set_ylim(0.5, 1.25)

    # Plot s2z_avg
    im2 = ax2_top.imshow(s2z_avg, cmap=balance, vmin=-1, vmax=1,
                        extent=extent, aspect="auto", origin="lower")
    # ax2_top.set_xlabel(L"\varphi")
    ax2_top.set_xticks([-pi/2, 0, pi/2])
    ax2_top.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax2_top.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax2_top.set_yticklabels([])
    ax2_top.set_ylim(0.5, 1.25)

    # Add shared colorbar for top row
    cbar_top = fig.colorbar(im1, ax=axs[1, :], orientation="vertical", pad=0.01)
    cbar_top.ax.xaxis.set_ticks_position("top")
    cbar_top.ax.xaxis.set_label_position("top")
    cbar_top.set_label(L"\langle s_z\rangle", 
                labelpad=0,  
                rotation=0,
                horizontalalignment="center",
                verticalalignment="bottom")
    cbar_top.ax.yaxis.label.set_position((0.25, 1.0))
    # Bottom row plots
    ax1_bottom, ax2_bottom = axs[2, 1], axs[2, 2]

    # Plot s1z_vrc
    im3 = ax1_bottom.imshow(s1z_vrc, cmap=lajolla.reversed(), vmin=0,
                            extent=extent, aspect="auto", origin="lower")
    ax1_bottom.set_ylabel(L"\eta/g")
    ax1_bottom.set_xticks([-pi/2, 0, pi/2])
    ax1_bottom.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax1_bottom.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax1_bottom.set_ylim(0.5, 1.25)

    # Plot s2z_vrc
    im4 = ax2_bottom.imshow(s2z_vrc, cmap=lajolla.reversed(), vmin=0,
                            extent=extent, aspect="auto", origin="lower")
    ax2_bottom.set_xlabel(L"\varphi")
    ax2_bottom.set_xticks([-pi/2, 0, pi/2])
    ax2_bottom.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax2_bottom.set_yticks([0.6, 0.8, 1.0, 1.2])
    ax2_bottom.set_yticklabels([])
    ax2_bottom.set_ylim(0.5, 1.25)

    # Add shared colorbar for bottom row
    cbar_bottom = fig.colorbar(im3, ax=axs[2, :], orientation="vertical", pad=0.01)
    cbar_bottom.ax.xaxis.set_ticks_position("top")
    cbar_bottom.ax.xaxis.set_label_position("top")
    cbar_bottom.set_label(L"\langle s_z^2\rangle", 
                labelpad=0,  
                rotation=0,
                horizontalalignment="center",
                verticalalignment="bottom")
    cbar_bottom.ax.yaxis.label.set_position((0.25, 1.0))

    return fig
end


function plot_spin_heatmaps(φ_vals, η_vals, s1x, s2x, s1y, s2y, s1z, s2z)
    # Create a 2x3 grid of subplots with a specified figure size
    fig, axs = subplots(2, 3, figsize=(15, 8))
    
    # Determine the extent for the heatmap axes (min and max of φ_vals and η_vals)
    extent = (minimum(φ_vals), maximum(φ_vals), minimum(η_vals), maximum(η_vals))
    
    # Plot Spin 1 components on the first row
    im1 = axs[1, 1].imshow(s1x, extent=extent, origin="lower", aspect="auto", cmap="viridis")
    axs[1, 1].set_title("Spin 1 x ")
    axs[1, 1].set_xlabel(L"\varphi")
    axs[1, 1].set_ylabel(L"\eta/g")
    colorbar(im1, ax=axs[1, 1])
    
    im2 = axs[1, 2].imshow(s1y, extent=extent, origin="lower", aspect="auto", cmap="viridis")
    axs[1, 2].set_title("Spin 1 y")
    axs[1, 2].set_xlabel(L"\varphi")
    axs[1, 2].set_ylabel(L"\eta/g")
    colorbar(im2, ax=axs[1, 2])
    
    im3 = axs[1, 3].imshow(s1z, extent=extent, origin="lower", aspect="auto", cmap="viridis")
    axs[1, 3].set_title("Spin 1 z")
    axs[1, 3].set_xlabel(L"\varphi")
    axs[1, 3].set_ylabel(L"\eta/g")
    colorbar(im3, ax=axs[1, 3])
    
    # Plot Spin 2 components on the second row
    im4 = axs[2, 1].imshow(s2x, extent=extent, origin="lower", aspect="auto", cmap="viridis")
    axs[2, 1].set_title("Spin 2 x")
    axs[2, 1].set_xlabel(L"\varphi")
    axs[2, 1].set_ylabel(L"\eta/g")
    colorbar(im4, ax=axs[2, 1])
    
    im5 = axs[2, 2].imshow(s2y, extent=extent, origin="lower", aspect="auto", cmap="viridis")
    axs[2, 2].set_title("Spin 2 y")
    axs[2, 2].set_xlabel(L"\varphi")
    axs[2, 2].set_ylabel(L"\eta/g")
    colorbar(im5, ax=axs[2, 2])
    
    im6 = axs[2, 3].imshow(s2z, extent=extent, origin="lower", aspect="auto", cmap="viridis")
    axs[2, 3].set_title("Spin 2 z")
    axs[2, 3].set_xlabel(L"\varphi")
    axs[2, 3].set_ylabel(L"\eta/g")
    colorbar(im6, ax=axs[2, 3])
    
    tight_layout()
    return fig
end

function compute_spin_observables(p::NamedTuple, adaga0::Int64, u0:: Vector{Float64}, tspan::Float64, w_cavity::Bool)
    K = p.κ / 2
    denominator = im * p.Δ - K
    # For two spins

    if w_cavity==true
        u0 = [u0..., sqrt(adaga0/2), sqrt(adaga0/2)]

        prob = ODEProblem(spin_cav_ode!, u0, tspan, p)

        sol = solve(prob, TRBDF2(autodiff=false), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)

    else
        prob = ODEProblem(spinonly_ode!, u0, tspan, p)

        sol = solve(prob, TRBDF2(autodiff=false), saveat=tspan/100,maxiters=1e9, abstol=1e-9, reltol=1e-9)
    end

    sx1_series = [u[1] for u in sol.u]
    sx2_series = [u[2] for u in sol.u]
    sy1_series = [u[3] for u in sol.u]
    sy2_series = [u[4] for u in sol.u]
    sz1_series = [u[5] for u in sol.u]
    sz2_series = [u[6] for u in sol.u]

    sx1_av, sx1_vrc = compute_time_average(sol.t, sx1_series)
    sx2_av, sx2_vrc = compute_time_average(sol.t, sx2_series)
    sy1_av, sy1_vrc = compute_time_average(sol.t, sy1_series)
    sy2_av, sy2_vrc = compute_time_average(sol.t, sy2_series)
    sz1_av, sz1_vrc = compute_time_average(sol.t, sz1_series)
    sz2_av, sz2_vrc = compute_time_average(sol.t, sz2_series)

    return  sx1_av, sx1_vrc, sx2_av, sx2_vrc, sy1_av, sy1_vrc, sy2_av, sy2_vrc, sz1_av, sz1_vrc, sz2_av, sz2_vrc
end

""" 
######################################################
FOURIER SPECTRUM 
######################################################
"""

# Function to compute FFT spectrum
function compute_fft_spectrum(signal, Δt)
    N = length(signal)
    Fs = 1 / Δt
    signal_centered = signal .- mean(signal)  # Remove DC component
    fft_vals = fft(signal_centered)
    freqs = fftfreq(N, Fs)
    # Shift frequencies and power to center DC component
    freqs_shifted = fftshift(freqs)
    power_shifted = fftshift(abs2.(fft_vals)) ./ (N * Fs)  # Two-sided PSD
    return freqs_shifted, power_shifted
end


#######PHASE DIAGRAM OF THE SPECTRUM########
function compute_dominant_freq(signal, Δt)
    freqs, power = compute_fft_spectrum(signal, Δt)
    Fs = 1/Δt
    nyquist = Fs / 2
    
    # Exclude DC and frequencies above Nyquist
    valid = findall(x -> 0 < abs(x) ≤ nyquist, freqs)
    power_valid = power[valid]
    freqs_valid = freqs[valid]
    
    # Thresholding: Ignore peaks below 10% of max power
    threshold = 0.1 * maximum(power_valid)
    significant = power_valid .> threshold
    if !any(significant)
        return NaN  # No meaningful peak
    end
    
    _, idx = findmax(power_valid[significant])
    return abs(freqs_valid[significant][idx])
end


function frequency_phase_diagram(φ_range, η_range, g, κ; show_progress=true)
    s1z_peaks = zeros(length(η_range), length(φ_range))
    s2z_peaks = similar(s1z_peaks)
    
    prog = Progress(length(φ_range) * length(η_range), enabled=show_progress)
    
    for (i, φ) in enumerate(φ_range)
        for (j, η) in enumerate(η_range)
            # Solve ODE (adjust parameters as per your system's requirements)
            p = (0.0, 0.0, 0.0, g, φ, κ, η)
            prob = ODEProblem(spinonly_ode!, [0.0, 0.0, 0.0, 0.0, -1.0, -1.0], (0.0, 10000.0), p)
            sol = solve(prob, Tsit5(), saveat=0.1, abstol=1e-8, reltol=1e-8)
            
            Δt = sol.t[2] - sol.t[1]
            s1z_peaks[j, i] = compute_dominant_freq(sol[5, :], Δt)
            s2z_peaks[j, i] = compute_dominant_freq(sol[6, :], Δt)
            next!(prog)
        end
    end
    return s1z_peaks, s2z_peaks
end 


function plot_frequency_phase_diagram(s1z_peaks, s2z_peaks, φ_range, η_range)
    zmin = min(minimum(s1z_peaks), minimum(s2z_peaks))
    zmax = max(maximum(s1z_peaks), maximum(s2z_peaks))
    clims = (zmin, zmax)


    # Create the figure and subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.4, 4), dpi=300)

    # Plot the first heatmap
    im1 = ax1.imshow(s1z_peaks, extent=[minimum(φ_range), maximum(φ_range), minimum(η_range ./ g), maximum(η_range ./ g)],
                    origin="lower", aspect="auto", cmap=lajolla.reversed(), vmin=clims[1], vmax=clims[2])
    ax1.set_xlabel(L"\varphi")
    ax1.set_ylabel(L"\eta/g")
    ax1.set_title(L"s_1")
    ax1.set_xticks([-π/2, 0, π/2])
    ax1.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax1.set_yticks([0.5, 1, 1.5])
    ax1.set_yticklabels(["0.5", "1.0", "1.5"])

    # Plot the second heatmap
    im2 = ax2.imshow(s1z_peaks, extent=[minimum(φ_range), maximum(φ_range), minimum(η_range ./ g), maximum(η_range ./ g)],
                    origin="lower", aspect="auto", cmap=lajolla.reversed(), vmin=zmin, vmax=zmax)
    ax2.set_xlabel(L"\varphi")
    ax2.set_title(L"s_2")
    ax2.set_xticks([-π/2, 0, π/2])
    ax2.set_xticklabels([L"-\pi/2", L"0", L"\pi/2"])
    ax2.set_yticks([0.5, 1, 1.5])
    ax2.set_yticklabels(["0.5", "1.0", "1.5"])

    # Create an axis for the colorbar
    cax = fig.add_axes([0.87, 0.22, 0.02, 0.65])  # [left, bottom, width, height]

    # Add the colorbar to the dedicated axis
    cbar = fig.colorbar(im1, cax=cax, orientation="vertical")
    cbar.set_label(L"\omega^\mathrm{max}", 
                labelpad=-20,  
                rotation=0,
                horizontalalignment="center",
                verticalalignment="bottom")
    # Manually position the label at the top using axes coordinates
    cbar.ax.yaxis.label.set_position((0.5, 1.0))
    plt.tight_layout(rect=[0, 0, 0.9, 1])

    return fig
end

function plot_time_and_freq_spectrum(t, signal1, freqs1, power1, signal2, freqs2, power2, g)
    fig = figure(figsize=(12, 8))

    # First subplot (top-left)
    subplot(2, 2, 1)
    plot(t, signal1, lw=5, color="steelblue")
    xlabel(L"\text{Time}")
    ylabel(L"s_1^z")  # Using LaTeX math mode
    
    
    # Second subplot (top-right)
    subplot(2, 2, 2)
    plot(freqs1, normalize(power1), lw=5, color="steelblue")
    xlabel(L"\text{Frequency}")
    ylabel(L"S(\omega)")
    xlim(-g/5, g/5)
    
    # Third subplot (bottom-left)
    subplot(2, 2, 3)
    plot(t, signal2, lw=5, color="coral")
    xlabel(L"\text{Time}")
    ylabel(L"s_2^z")
    
    
    # Fourth subplot (bottom-right)
    subplot(2, 2, 4)
    plot(freqs2, normalize(power2), lw=5, color="coral")
    xlabel(L"\text{Frequency}")
    ylabel(L"S(\omega)")
    xlim(-g/5, g/5)
    
    tight_layout()
    return fig
end

""" 
############################
RADIANCE WITNESS 
############################
"""

function compute_cavity_observables(p::NamedTuple, adaga0::Int64, u0_two:: Vector{Float64}, u0_one:: Vector{Float64}, tspan::Float64, w_cavity::Bool)
    K = p.κ / 2
    denominator = im * p.Δ - K
    # For two spins

    if w_cavity==true
        u0_two = [u0_two..., sqrt(adaga0/2), sqrt(adaga0/2)]
        u0_one = [u0_one..., sqrt(adaga0/2), sqrt(adaga0/2)]
        prob_two = ODEProblem(spin_cav_ode!, u0_two, tspan, p)
        prob_one = ODEProblem(single_spin_cav_ode!, u0_one, tspan, p)
        sol_two = solve(prob_two, TRBDF2(autodiff=false), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)
        sol_one = solve(prob_one, TRBDF2(autodiff=false), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)

        # For two spins
        a_two = sol_two[7, :] .+ im * sol_two[8, :]
        adaga_two = abs2.(a_two)
        avg_two = compute_time_average(sol_two.t, adaga_two)

        # For one spin
        a_one = sol_one[4, :] .+ im * sol_one[5, :]
        adaga_one = abs2.(a_one)
        avg_one = compute_time_average(sol_one.t, adaga_one)
    else
        prob_two = ODEProblem(spinonly_ode!, u0_two, tspan, p)
        prob_one = ODEProblem(single_spinonly_ode!, u0_one, tspan, p)
        sol_two = solve(prob_two, TRBDF2(autodiff=false), saveat=tspan/100,maxiters=1e9, abstol=1e-9, reltol=1e-9)
        sol_one = solve(prob_one, TRBDF2(autodiff=false), saveat=tspan/100, maxiters=1e9, abstol=1e-9, reltol=1e-9)
        
        # For two spins
        s1x_two = sol_two[1, :]
        s2x_two = sol_two[2, :]
        s1y_two = sol_two[3, :]
        s2y_two = sol_two[4, :]
    
        a_two = @. ( ( g/2 * (im*(s1x_two + exp(-im*φ)*s2x_two) + (s1y_two + exp(-im*φ)*s2y_two)) - η ) / denominator )
        adaga_two = abs2.(a_two)
        avg_two = compute_time_average(sol_two.t, adaga_two)
    
        # For one spin    
        s1x_one = sol_one[1, :]
        s1y_one = sol_one[2, :]
    
        a_one = @. ( (g/2 * ((im*s1x_one + s1y_one)) - η ) / denominator )
        adaga_one = abs2.(a_one)
        avg_one = compute_time_average(sol_one.t, adaga_one)
    end
    # Compute R
    R = (avg_two -  avg_one) / ( avg_one)
    return R, avg_two, avg_one
end


"""
############################################
MAIN SCRIPT (NUMERICS)
############################################    
"""
φ = 0
η = κ
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)

t_end = 100/κ
u0 = [0.0, 0.0,0.0, 0.0,-1.0, -1.0, 0, 0]

prob = ODEProblem(spin_cav_ode!, u0, (0.0, t_end), p)
sol = solve(prob, TRBDF2(autodiff=false); saveat=t_end/100)

sx1 = 
s1sq = sol[1].^2+s1y^2+s1z^2

######################################################
# Parameters
Δ1 = 0
Δ2 = 0
Δ = 0
g = 1.0
κ = 1.0
φ_range = range(-π/2, π/2, 41);    # φ (x-axis)
η_range = range(0.5*g, 1.5*g, 41);  # η (y-axis)

u0 = [0.0, 0.0,    # spin 1 (e.g. sx, sy)
      0.0, 0.0,    # spin 2 (e.g. sx, sy)
      -1.0, -1.0];    # z-components
w_cavity = true ;
adaga0=0;
tspan=100/κ;

s1x = zeros(length(η_range), length(φ_range));
s2x = zeros(length(η_range), length(φ_range));
s1y = zeros(length(η_range), length(φ_range));
s2y = zeros(length(η_range), length(φ_range));
s1z = zeros(length(η_range), length(φ_range));
s2z = zeros(length(η_range), length(φ_range));

# ---------------------------------------------------
# Set up progress bar over the total number of simulations
# ---------------------------------------------------
n_total = length(η_range) * length(φ_range);
progress = Progress(n_total, desc="Evaluating...")

# ---------------------------------------------------
# Loop over the grid of φ and η values
# ---------------------------------------------------
for (i, η_val) in enumerate(η_range)
    for (j, φ_val) in enumerate(φ_range)
        # Create the parameter tuple with the current φ and η values
        p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ_val, κ=κ, η=η_val)
        
        # Use a fresh copy of u0 to avoid modifying the original
        sx1_av, sx1_vrc, sx2_av, sx2_vrc, sy1_av, sy1_vrc, sy2_av, sy2_vrc, sz1_av, sz1_vrc, sz2_av, sz2_vrc =
            compute_spin_observables(p, adaga0, u0, tspan, w_cavity)
        
        # Store the average values for each observable
        s1x[i, j] = sx1_av
        s2x[i, j] = sx2_av
        s1y[i, j] = sy1_av
        s2y[i, j] = sy2_av
        s1z[i, j] = sz1_av
        s2z[i, j] = sz2_av
        
        # Update the progress bar
        next!(progress)
    end
end
fig, axs = subplots(2, 3, figsize=(15, 10))

function plot_heatmap(ax, data, title_str)
    im = ax.imshow(data, extent=(first(φ_range), last(φ_range), first(η_range), last(η_range)),
                   origin="lower", aspect="auto")
    ax.set_title(title_str)
    ax.set_xlabel(L"\varphi")
    ax.set_ylabel(L"\eta")
    fig.colorbar(im, ax=ax)
end

plot_heatmap(axs[1, 1], s1x, L"s1x")
plot_heatmap(axs[2, 1], s2x,  L"s2x")
plot_heatmap(axs[1, 2], s1y,  L"s1y")
plot_heatmap(axs[2, 2], s2y,  L"s2y")
plot_heatmap(axs[1, 3], s1z,  L"s1z")
plot_heatmap(axs[2, 3], s2z, L"s2z")

tight_layout()
display(gcf())



spin1_length = s1x.^2 .+ s1y.^2 .+ s1z.^2;
spin2_length = s2x.^2 .+ s2y.^2 .+ s2z.^2;

# -----------------------------
# Plot the spin lengths in a new figure with two subplots
# -----------------------------
fig2, axs2 = subplots(1, 2, figsize=(12, 5));

im1 = axs2[1].imshow(spin1_length, extent=(first(φ_range), last(φ_range), first(η_range), last(η_range)),
                      origin="lower", aspect="auto");
axs2[1].set_title(L"s_1^2");
axs2[1].set_xlabel(L"\varphi");
axs2[1].set_ylabel(L"\eta");
fig2.colorbar(im1, ax=axs2[1]);

im2 = axs2[2].imshow(spin2_length, extent=(first(φ_range), last(φ_range), first(η_range), last(η_range)),
                      origin="lower", aspect="auto");
axs2[2].set_title(L"s_2^2")
axs2[2].set_xlabel(L"\varphi");
axs2[2].set_ylabel(L"\eta");
fig2.colorbar(im2, ax=axs2[2]);

tight_layout()
display(gcf())

###################################################################

Δ1 = 0
Δ2 = 0
Δ = 0 
κ = 1.0
g = 1κ
φ = 0
η = κ
t_end = 10/κ
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)
saveat = t_end/1000

u0c_coupled=[0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0, 0]
u0s_coupled=[0.0, 0.0, 0.0, 0.0, -1.0, -1.0]
u0c_single=[0.0, 0.0, -1.0, 0, 0]
u0s_single=[0.0, 0.0, -1.0]

function spin_cav_dyns(p::NamedTuple,t_end, saveat, u0c_coupled,u0s_coupled,u0c_single,u0s_single,
    maxiters=1e9, reltol=1e-9, abstol=1e-9)

    # Precompute terms used in post-processing
    K = p.κ / 2
    denominator = im * p.Δ - K
    pre_g = (p.g / 2) / denominator
    pre_η = p.η / denominator
    φ_term = exp(-im * p.φ)

    # Common solver options
    solver_opts = (maxiters=maxiters, reltol=reltol, abstol=abstol)

    # 1. Coupled Spins WITH Cavity
    prob_cav_coupled = ODEProblem(spin_cav_ode!, u0c_coupled, (0.0, t_end), p)
    sol_cav_coupled = solve(prob_cav_coupled, Tsit5(); saveat=saveat, solver_opts...)

    # 2. Coupled Spins WITHOUT Cavity
    prob_spin_coupled = ODEProblem(spinonly_ode!, u0s_coupled, (0.0, t_end), p)
    sol_spin_coupled = solve(prob_spin_coupled, Tsit5(); saveat=saveat, solver_opts...)

    # 3. Single Spin WITH Cavity
    prob_cav_single = ODEProblem(single_spin_cav_ode!, u0c_single, (0.0, t_end), p)
    sol_cav_single = solve(prob_cav_single, Tsit5(); saveat=saveat, solver_opts...)

    # 4. Single Spin WITHOUT Cavity
    prob_spin_single = ODEProblem(single_spinonly_ode!, u0s_single, (0.0, t_end), p)
    sol_spin_single = solve(prob_spin_single, Tsit5(); saveat=saveat, solver_opts...)

    # Calculate |a|² for all cases
    abs2_cav_coupled = abs2.(sol_cav_coupled[7,:] .+ im.*sol_cav_coupled[8,:])

    s1x, s2x = sol_spin_coupled[1,:], sol_spin_coupled[2,:]
    s1y, s2y = sol_spin_coupled[3,:], sol_spin_coupled[4,:]
    abs2_spin_coupled = abs2.(pre_g .* (im.*(s1x .+ φ_term.*s2x) .+ (s1y .+ φ_term.*s2y)) .- pre_η)

    abs2_cav_single = abs2.(sol_cav_single[4,:] .+ im.*sol_cav_single[5,:])

    s1x, s1y = sol_spin_single[1,:], sol_spin_single[2,:]
    abs2_spin_single = abs2.(pre_g .* (im.*s1x .+ s1y) .- pre_η)

    avg_twoc = compute_time_average(sol_spin_coupled.t, abs2_cav_coupled)
    avg_onec = compute_time_average(sol_spin_single.t, abs2_cav_single)
    R_cav = ( avg_twoc- 2 * avg_onec) / (2 * avg_onec)

    avg_twos = compute_time_average(sol_spin_coupled.t, abs2_spin_coupled)
    avg_ones = compute_time_average(sol_spin_single.t, abs2_spin_single)
    R_spin = ( avg_twos- 2 * avg_ones) / (2 * avg_ones)


    R_cav_str = string(round(R_cav, digits=3))
    R_spin_str = string(round(R_spin, digits=3))
    # Create plot
    fig = figure(figsize=(8,5))
    t_points = sol_cav_coupled.t

    plot(t_points, abs2_cav_coupled, lw=2.5, color="dodgerblue", label="Coupled spins (with cavity)")
    #plot(t_points, abs2_spin_coupled, lw=2.5, color="dodgerblue", linestyle="--", label="Coupled spins (spin-only)")
    plot(t_points, abs2_cav_single, lw=2.5, color="crimson", label="Single spin (with cavity)")
    #plot(t_points, abs2_spin_single, lw=2.5, color="crimson", linestyle="--", label="Single spin (spin-only)")

    anon1 = annotate(L"$\mathcal{R}_c = " * R_cav_str * L"$", 
    (0.5, 0.7), 
    xycoords="axes fraction",
    fontsize=18,
    color="black")

    anon2 = annotate(L"$\mathcal{R}_s = " * R_spin_str * L"$", 
    (0.5, 0.6), 
    xycoords="axes fraction",
    fontsize=18,
    color="black")

    xlabel(L"Time (1/$\kappa$)")
    ylabel(L"|\alpha|^2")
    legend()
    grid(true)
    tight_layout()

    return fig
end


compute_time_average(sol_cav_coupled.t, abs2_cav_coupled)
compute_time_average(sol_spin_coupled.t, abs2_spin_coupled)

avg_two = compute_time_average(sol_spin_coupled.t, abs2_cav_coupled)
avg_one = compute_time_average(sol_spin_single.t,  abs2_cav_single)

R = (avg_two - 2 * avg_one) / (2 * avg_one)


# Extract time series (using second half)
t = sol.t
signal1 = sol[5, :]
signal2 = sol[6, :]



# Compute spectra
Δt = t[2] - t[1]
freqs1, power1 = compute_fft_spectrum(signal1, Δt)
freqs2, power2 = compute_fft_spectrum(signal2, Δt)

fig = plot_time_and_freq_spectrum(t, signal1, freqs1, power1, signal2, freqs2, power2, g)

###################################################################


"""
############
 eta vs g phase diag
############
 """
κ = 1
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)

η_range = range(0.5*κ, 5κ, length=21)  
g_range = range(0.5*κ, 5κ, length=21)

# Define the range of φ values
φ_values = [0, π/4, π/2, π];

# Initialize grids for each φ value
adaga2_grids = [];
adaga1_grids = [];
# Total iterations for the φ loop
total_phi_iterations = length(φ_values);
prog_phi = Progress(total_phi_iterations, 1, "Computing φ grids: ", 50);

# Loop over φ values and compute the phase diagrams
for φ in φ_values
    # Update the parameter set with the current φ
    p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)

    # Initialize grids
    adaga2_grid = zeros(length(η_range), length(g_range))
    adaga1_grid = zeros(length(η_range), length(g_range))
    # Progress bar for the η, g loop
    total_iterations = length(η_range) * length(g_range)
    prog = Progress(total_iterations, 1, "Computing η, g grid for φ = $φ: ", 50)

    # Compute the phase diagram for the current φ
    for (i, η) in enumerate(η_range)
        for (j, g) in enumerate(g_range)
            # Create new parameter set with current g and η
            current_p = merge(p, (g=g, η=η))
            # Compute cavity observables
            _, adaga2, adaga1 = compute_cavity_observables(current_p, 0, u0_two, u0_one, 100.0/κ, true)
            adaga2_grid[i, j] = adaga2
            adaga1_grid[i, j] = adaga1
            next!(prog)  # Update the η, g progress bar
        end
    end

    # Store the grid for the current φ
    push!(adaga2_grids, adaga2_grid)
    push!(adaga1_grids, adaga1_grid)
    next!(prog_phi)  # Update the φ progress bar
end

# Plot the 2x2 grid of phase diagrams
fig, axs = subplots(2, 2, figsize=(12, 10));

# Titles for each subplot
titles = [L"\varphi = 0", L"\varphi = \pi/4", L"\varphi = \pi/2", L"\varphi = \pi"];

img = nothing
# Plot each phase diagram
for (idx, ax) in enumerate(axs)
    φ_idx = div(idx - 1, 2) + 1
    η_idx = mod(idx - 1, 2) + 1
    grid_idx = (φ_idx - 1) * 2 + η_idx
    img = ax.imshow(adaga2_grids[grid_idx], extent=[minimum(g_range), maximum(g_range), minimum(η_range), maximum(η_range)],
              origin="lower", aspect="auto", cmap=cet_cbl, norm=LogNorm(vmin=1e-3, vmax=maximum(adaga2_grids[grid_idx])))
    ax.set_title(titles[grid_idx])
    if idx == 1 || idx==3 
        ax.set_xlabel("")
    else
        ax.set_xlabel(L"g/\kappa")  # Only set for the bottom row
    end

    # Remove y-axis labels for the right-hand column
    if idx == 3 || idx== 4 
        ax.set_ylabel("")
    else
        ax.set_ylabel(L"\eta/\kappa")  # Only set for the left column
    end
end

fig.subplots_adjust(right=0.85);
cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.7]) ;
cbar = fig.colorbar(img, cax=cbar_ax);
cbar.set_label(L"|\alpha|^2");
# Adjust layout and save the figure
#tight_layout()
savefig("eta_g_adaga_zoom_t100.pdf", dpi=300)
display(fig)




g = 1
κ = 1
η = 0
φ = pi
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)



η_range = range(0.1*g, 2*g, length=1000);  
total_iterations = length(η_range);
prog = Progress(total_iterations, 1, "Computing: ", 50);
adaga2_grid = zeros(length(η_range));
adaga1_grid  = zeros(length(η_range));
R_grid = zeros(length(η_range));
for (j, η) in enumerate(η_range)
    current_p = merge(p, (η=η,)) # <- Note the comma to create a NamedTuple
    R, adaga2, adaga1 = compute_cavity_observables(current_p, 0, u0_two, u0_one, 10.0/κ, true)
    adaga2_grid[j] = adaga2
    adaga1_grid[j] = adaga1
    R_grid[j] = R
    next!(prog)
end

fig = figure(figsize=(8, 6))
ax = gca()

# Plot |α|² vs η/g
fig, ax = plt.subplots(figsize=(8, 6))
plot(η_range ./ g, adaga2_grid, "b-", linewidth=2, label="Coupled")
plot(η_range ./ g, adaga1_grid, "r-", linewidth=2, label="Single")
plot(η_range ./ g, R_grid, "g-", linewidth=2, label="Radiance")
xlabel(L"\eta / g")
ylabel(L"|\alpha|^2")
grid(true)
tight_layout()
legend()
display(fig)



g = 1.0
κ = 1
η = 0
φ =  0
p = (Δ1=Δ1, Δ2=Δ2, Δ=Δ, g=g, φ=φ, κ=κ, η=η)

η_range = range(0.4g, 1.5g, length=51)  
ϕ_range = range(-pi/2, pi/2, length=51)
total_iterations = length(η_range) * length(ϕ_range)
prog = Progress(total_iterations, 1, "Computing grids: ", 50);
R_grid = zeros(length(η_range), length(ϕ_range));
adaga2_grid = zeros(length(η_range), length(ϕ_range));
adaga1_grid = zeros(length(η_range), length(ϕ_range));
for (i, η) in enumerate(η_range)
    for (j, φ) in enumerate(ϕ_range)
        # Create new parameter set with current phi and eta
        current_p = merge(p,(η=η, φ=φ))
        # Compute R for this parameter set
        R, adaga2, adaga1 = compute_cavity_observables(current_p, 0, u0_two, u0_one, 100.0/κ, true)
        R_grid[i, j] = R
        adaga2_grid[i, j] = adaga2
        adaga1_grid[i, j] = adaga1
        next!(prog)
    end
end


fig = figure(figsize=(8, 6));   
ax = gca();
extent = [minimum(ϕ_range)/π, maximum(ϕ_range)/π,minimum(η_range)/g, maximum(η_range)/g];
img_plot = ax.imshow(adaga2_grid, extent=extent,origin="lower",aspect="auto",cmap=balance.reversed(), norm=LogNorm(vmin=1e-3, vmax=1e2));
cbar = colorbar(img_plot, ax=ax);
#cbar.set_label(L"\Delta^\alpha_{1,2}");
cbar.set_label(L"|\alpha|^2");

#cbar.set_label(L"\mathcal{R}");
xlabel(L"\varphi/\pi")
ylabel(L"\eta/g")
tight_layout()
display(fig)
