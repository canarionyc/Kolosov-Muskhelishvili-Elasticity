# %%
# import Pkg; Pkg.add("Symbolics")
# import Pkg; Pkg.add("LinearAlgebra")
using Symbolics, LinearAlgebra

# ==============================================================================
## 1. SETUP: Define Variables and Coordinates
# ==============================================================================
@variables r theta phi          # Coordinates
@variables U(r)                 # Radial displacement function
@variables λ μ                  # Lame's Elastic Constants (Love uses λ, μ)

# Define the Metric Tensor g_ij for Spherical Coordinates (r, θ, φ)
# ds^2 = 1*dr^2 + r^2*dθ^2 + r^2*sin(θ)^2*dφ^2
g_matrix = Diagonal([1, r^2, r^2 * (sin(theta))^2])
g_inv = inv(g_matrix)

# ==============================================================================
## 2. PHYSICS ENGINE: Automated Tensor Calculus
# ==============================================================================
# Function to compute Christoffel Symbols of the Second Kind (Γ^k_ij)
# Γ^k_ij = 0.5 * g^kl * (∂g_lj/∂x_i + ∂g_li/∂x_j - ∂g_ij/∂x_l)
function compute_christoffel(g, g_inv, vars)
    n = length(vars)
    Γ = zeros(Num, n, n, n)
    for k in 1:n, i in 1:n, j in 1:n
        sum_term = 0
        for l in 1:n
            # Partial derivatives of the metric
            dg_lj_i = Symbolics.derivative(g[l, j], vars[i])
            dg_li_j = Symbolics.derivative(g[l, i], vars[j])
            dg_ij_l = Symbolics.derivative(g[i, j], vars[l])
            sum_term += g_inv[k, l] * (dg_lj_i + dg_li_j - dg_ij_l)
        end
        Γ[k, i, j] = 0.5 * sum_term
    end
    return Γ
end

vars = [r, theta, phi]
Γ = compute_christoffel(g_matrix, g_inv, vars)

# ==============================================================================
## 3. KINEMATICS: Displacement and Strain
# ==============================================================================
# Define Displacement Vector u^i = [U(r), 0, 0] (Purely radial)
u = [U(r), 0, 0]

# Compute Linearized Strain Tensor: ε_ij = 1/2 (∇_i u_j + ∇_j u_i)
# We lower indices first: u_j = g_jk * u^k
u_cov = g_matrix * u

# Covariant Derivative: ∇_i u_j = ∂u_j/∂x_i - Γ^k_ij * u_k
n = 3
grad_u = Matrix{Num}(undef, n, n)
for i in 1:n, j in 1:n
    term1 = Symbolics.derivative(u_cov[j], vars[i])
    term2 = sum(Γ[k, i, j] * u_cov[k] for k in 1:n)
    grad_u[i, j] = term1 - term2
end

strain_tensor = 0.5 .* (grad_u + grad_u')

# ==============================================================================
# 4. CONSTITUTIVE LAW: Stress (Isotropic)
# ==============================================================================
# Trace of Strain (Dilatation Δ)
# Δ = g^ij * ε_ij
dilatation = sum(g_inv[i,j] * strain_tensor[i,j] for i in 1:n, j in 1:n)

# Stress Tensor σ^ij = λ * Δ * g^ij + 2μ * g^ik * g^jl * ε_kl
# (Using contravariant form for easier divergence calculation later)
stress_contra = Matrix{Num}(undef, n, n)
for i in 1:n, j in 1:n
    term_lambda = λ * dilatation * g_inv[i, j]
    term_mu = 0
    for k in 1:n, l in 1:n
        term_mu += 2 * μ * g_inv[i, k] * g_inv[j, l] * strain_tensor[k, l]
    end
    stress_contra[i, j] = term_lambda + term_mu
end

# ==============================================================================
# 5. EQUILIBRIUM: Divergence of Stress = 0
# ==============================================================================
# Only the radial component (i=1) is non-trivial due to symmetry
# ∇_j σ^1j = ∂σ^1j/∂x_j + Γ^1_jk σ^jk + Γ^j_jk σ^1k
div_sigma_r = 0
i = 1 # Radial index
for j in 1:n
    # Partial derivative
    d_sigma = Symbolics.derivative(stress_contra[i, j], vars[j])
    
    # Christoffel terms
    gamma_term1 = sum(Γ[i, j, k] * stress_contra[j, k] for k in 1:n)
    gamma_term2 = sum(Γ[j, j, k] * stress_contra[i, k] for k in 1:n)
    
    global div_sigma_r += d_sigma + gamma_term1 + gamma_term2
end

# Simplify the result
equilibrium_eq = simplify(div_sigma_r)

println("Calculated Equilibrium Equation (Radial):")
println(equilibrium_eq)