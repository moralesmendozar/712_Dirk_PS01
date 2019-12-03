# -----------------------------------------------------------------------------
# DISCRETIZATION OF A CONTINUOUS STOCHASTIC PROCESS
# Methods:
# (a) Tauchen Method:
# George Tauchen (1986): "Finite state markov-chain approximations to univariate and vector autoregressions", Economic Letters
#
# (b) Rouwenhorst Method:
# Geert Rouwenhorst (1995): "Asset pricing implications of equilibirium business cycle models", in Thomas F.Cooley "Frontiers of business cycle research"
#
# Coded by Juan Castellanos Silván using Spencer Lyon's rutines in QuantEcon
# October 20th, 2019
# Julia, Version 1.2.0
# ------------------------------------------------------------------------------
module discretizingMethods

export tauchen, rouwenhorst

    using Optim
    using Distributions: cdf, Normal, quantile

    """
    Tauchen's (1996) method for approximating AR(1) process with finite markov chain. The process follows:
        y_t = (1-ρ)μ + ρ y_{t-1} + ε_t,   where ε_t ∼ N (0, σ^2)

    Recall that E[y_t] = μ and Var[y_t] = σ_y = σ^2 / (1-ρ^2)

    # ----------
    # Arguments
    # ----------
    - `N::Integer`: Number of points in markov process
    - `ρ::Real` : Persistence parameter in AR(1) process
    - `σ::Real` : Standard deviation of random component of AR(1) process
    - `μ::Real(0.0)` : Mean of AR(1) process
    - `m::Integer(3)` : The number of standard deviations to each side the process should span

    # --------
    # Returns
    # -------
    - `mTranstion::Array{Float64,2}` : N x N Markov transition matrix
    - `mStationary::Array{Float64,2}`: N x 1 Stationary Markov chain
    - `vGrid::Array{Float64,1}`: N vector with discrete state values
    """
    function tauchen(N::Integer, ρ::Real, σ::Real, μ::Real=0.0, m::Integer=3)

        # Discretize state space
        σ_y = σ / sqrt(1-ρ^2) # standard deviation of y_t
        y_max = μ + m * σ_y # upper bound
        y_min = μ - m * σ_y # lower bound
        d = (y_max - y_min)/(N-1) # length of the interval
        vGrid = collect(y_min:d:y_max)

        # Transition matrix
        mTransition = zeros(N,N)
        for iRow in 1:N
            for jCol in 2:N-1
                mTransition[iRow,jCol]=cdf(Normal(), (vGrid[jCol] + d/2 - (1-ρ)*μ - ρ*vGrid[iRow]) / σ) - cdf(Normal(), (vGrid[jCol] - d/2 - (1-ρ)*μ - ρ*vGrid[iRow]) / σ)
            end
            mTransition[iRow,1] = cdf(Normal(), (vGrid[1] + d/2 - (1-ρ)*μ - ρ*vGrid[iRow]) / σ)
            mTransition[iRow,N] = 1 - cdf(Normal(),(vGrid[N] - d/2 - (1-ρ)*μ - ρ*vGrid[iRow]) / σ)
        end

        if !(sum(mTransition, dims = 2) ≈ ones(N,1))
            error("At least one of the rows do not sum up to 1")
        end

        return vGrid, mTransition
    end


    """
    Rouwenhorst's method to approximate AR(1) processes.
    The process follows:
        y_t = ρ y_{t-1} + ε_t,  where ε_t ∼ N (0, σ^2)

    # ----------
    # Arguments
    # ----------
    - `N::Integer` : Number of points in markov process
    - `ρ::Real` : Persistence parameter in AR(1) process
    - `σ::Real` : Standard deviation of random component of AR(1) process
    - `μ::Real(0.0)` :  Mean of AR(1) process

    # --------
    # Returns
    # -------
    - `mTranstion::Array{Float64,2}` : N x N Markov transition matrix
    - `mStationary::Array{Float64,2}`: N x 1 Stationary Markov chain
    - `vGrid::Array{Float64,1}`: N vector with discrete state values
    """
    function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)

        σ_y = σ / sqrt(1-ρ^2)
        p  = (1+ρ)/2
        Θ = [p 1-p; 1-p p]
        ψ = sqrt(N-1) * σ_y
        m = μ / (1 - ρ)

        vGrid, mTransition = _rouwenhorst(p,p,m,ψ,N)

        return collect(vGrid), mTransition
    end

    function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
        if n == 2
            return [m-Δ, m+Δ], [p 1-p; 1-q q]
        else
            _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n-1)

            θN = p* [θ_nm1 zeros(n-1, 1); zeros(1, n)] +
                (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
                q*[zeros(1, n); zeros(n-1, 1) θ_nm1] +
                (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]

            θN[2:end-1, :] ./= 2

            return range(m-Δ, stop=m+Δ, length=n), θN
        end
    end
end
