"""
    waic(ll::AbstractArray{<:Real}; pointwise=false, log_lik="log_lik, kwargs...)

Compute the Widely Applicable Information Criterion (WAIC).

# Arguments
* `ll::AbstractMatrix`      : A matrix of posterior log likelihoods with (samples, observations) shape
* `pointwise::Bool`         : Compute WAIC pointwise, return a vector


# Returns
* `res::NamedTuple`: (WAIC=waics, lppd=lpd, penalty=pD, std_err=se) where

    WAIC                    : Sum of pointwise waic values (or pointwise vector)
    lppd                    : Log pointwise predictive density
    penalty                 : Penalty term ("overfitting penalty")
    std_err                 : Standard error of pointwise waic values

"""
function waic(ll::AbstractMatrix; pointwise=false)
    n_samples, n_obs = size(ll)
    pD = var2.(eachcol(ll))
    lpd = reshape(logsumexp(ll .- log(n_samples); dims=1), n_obs);

    waic_vec = -2 * (lpd - pD);
    if pointwise
        waics = waic_vec
    else
        waics = sum(waic_vec)
        lpd = sum(lpd)
        pD = sum(pD)
    end

    se = n_obs*var2(waic_vec)
    se = se < 0.0 ? nothing : sqrt(se)
    (WAIC=waics, lppd=lpd, penalty=pD, std_err=se)
end

export waic
