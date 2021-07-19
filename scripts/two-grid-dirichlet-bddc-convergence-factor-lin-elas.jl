# dependencies
using LFAToolkit
using LinearAlgebra
using CSV
using DataFrames

# setup
dimension = 3
numbercomponents = 3
mesh = Mesh3D(1.0, 1.0, 1.0)
convergencefactors = DataFrame()

# test range
for finep = 2:4
    println("fine_p = ", finep)
    for coarsep = 1:finep-1
        println("  coarse_p = ", coarsep)
        # setup
        # -- bases
        finebasis = TensorH1LagrangeBasis(finep + 1, finep + 1, numbercomponents, dimension)
        coarsebasis =
            TensorH1LagrangeBasis(coarsep + 1, finep + 1, numbercomponents, dimension)
        ctofbasis = TensorH1LagrangeBasis(
            coarsep + 1,
            finep + 1,
            numbercomponents,
            dimension,
            lagrangequadrature = true,
        )

        # constants
        e = 1E6                             # Young's modulus
        ν = 0.3                             # Poisson's ratio
        λ = e * ν / ((1 + ν) * (1 - 2 * ν)) # Lamé parameters
        μ = e / (2 * (1 + ν))

        function linearelasticityweakform(deltadu::Array{Float64}, w::Array{Float64})
            # strain
            dϵ = (deltadu + deltadu') / 2
            # strain energy
            dσ_11 = (λ + 2μ) * dϵ[1, 1] + λ * dϵ[2, 2] + λ * dϵ[3, 3]
            dσ_22 = λ * dϵ[1, 1] + (λ + 2μ) * dϵ[2, 2] + λ * dϵ[3, 3]
            dσ_33 = λ * dϵ[1, 1] + λ * dϵ[2, 2] + (λ + 2μ) * dϵ[3, 3]
            dσ_12 = μ * dϵ[1, 2]
            dσ_13 = μ * dϵ[1, 3]
            dσ_23 = μ * dϵ[2, 3]
            dσ = [dσ_11 dσ_12 dσ_13; dσ_12 dσ_22 dσ_23; dσ_13 dσ_23 dσ_33]

            # delta dv
            deltadv = dσ * w[1]

            return [deltadv']
        end

        # linear elasticity operators
        function makeoperator(basis::TensorBasis)
            inputs = [
                OperatorField(basis, [EvaluationMode.gradient], "gradent of deformation"),
                OperatorField(
                    basis,
                    [EvaluationMode.quadratureweights],
                    "quadrature weights",
                ),
            ]
            outputs = [
                OperatorField(
                    basis,
                    [EvaluationMode.gradient],
                    "test function gradient of deformation",
                ),
            ]
            return Operator(linearelasticityweakform, mesh, inputs, outputs)
        end
        fineoperator = makeoperator(finebasis)
        coarseoperator = makeoperator(coarsebasis)

        # -- smoothers
        identity = IdentityPC(fineoperator)
        bddc = DirichletBDDC(fineoperator)

        # -- p-multigrid preconditioner
        multigrid = PMultigrid(fineoperator, coarseoperator, identity, [ctofbasis])

        # compute smoothing factor
        # -- setup
        numbersteps = 8
        maxeigenvalue = 0
        θ_min = -π / 2
        θ_min_high = π / 2
        θ_max = 3π / 2
        θ_step = 2π / numbersteps
        θ_range = θ_min:θ_step:(θ_max-θ_step)
        ω_range = 0.01:0.01:0.29
        num_ω = max(size(ω_range)...)
        v_range = 1:1
        num_v = 1

        # -- compute
        maxeigenvalue = zeros(num_ω, num_v)
        θ_maxeigenvalue = -1 * ones(num_ω, num_v, dimension)
        # -- 3D --
        for i = 1:numbersteps, j = 1:numbersteps, k = 1:numbersteps
            θ = [θ_range[i], θ_range[j], θ_range[k]]
            if sqrt(abs(θ[1])^2 + abs(θ[2])^2 + abs(θ[3])^2) > π / 128
                M = computesymbols(multigrid, [0], [0, 0], θ)
                A = I - computesymbols(bddc, [1.0], θ)
                for w = 1:num_ω, v in v_range
                    S = (I - ω_range[w] * A)^v
                    eigenvalues = [abs(val) for val in eigvals(S * M * S)]
                    currentmaxeigenvalue = max(eigenvalues...)
                    if (currentmaxeigenvalue > maxeigenvalue[w, v])
                        maxeigenvalue[w, v] = currentmaxeigenvalue
                        θ_maxeigenvalue[w, v, :] = θ / π
                    end
                end
            end
        end
        for v in v_range
            (_, ω_index) = findmin(maxeigenvalue[:, v])
            append!(
                convergencefactors,
                DataFrame(
                    finep = finep,
                    coarsep = coarsep,
                    ω = ω_range[ω_index],
                    θ = θ_maxeigenvalue[ω_index, v, :],
                    v = v,
                    ρ = maxeigenvalue[ω_index, v],
                ),
            )
        end
    end
end

CSV.write("two-grid-dirichlet-bddc-convergence-factor-lin-elas.csv", convergencefactors)
