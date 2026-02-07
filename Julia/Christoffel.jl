## 
import Pkg; Pkg.add("ForwardDiff")

using ForwardDiff  # For automatic differentiation of the metric

## Define a coordinate transformation: (r, θ, φ) -> (x, y, z)
# Example: Spherical Coordinates
function transformation(q)
    r, θ, ϕ = q
    return [r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ)]
end

# Compute the Metric Tensor g_ij = J^T * J
function get_metric(q)
    J = ForwardDiff.jacobian(transformation, q)
    return J' * J
end

## Compute Christoffel Symbols (Γ) using the metric
# Γ^k_ij = 0.5 * g^kl * (∂g_lj/∂q_i + ∂g_li/∂q_j - ∂g_ij/∂q_l)
# This replaces Love's manual "scale factor" geometry.

@variables r θ ϕ
get_metric([r, θ, ϕ])