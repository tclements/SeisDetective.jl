export recursive_sta_lta, classic_sta_lta, carl_sta_trig, delayed_sta_lta
export z_detect, trigger_onset

"""
  recursive_sta_lta(A, nsta, nlta)

Computes recursive STA/LTA. 

The length of the STA and LTA are given by `nsta` and `lsta` in samples, respectively.

# Arguments
- `A::AbstractArray`: 1D array.
- `nsta:Int`: Length of short time average window in samples
- `nlta:Int`: Length of long time average window in samples
"""
function recursive_sta_lta(A::AbstractArray, nsta::Int, nlta::Int)
    N = size(A,1)

    # compute short time average (STA) and long term average (LTA)
    # given by Evans and Allen 
    csta = 1.0 / nsta 
    clta = 1.0 / nlta 
    sta = 0.0 
    lta = 1e-99 # avoid zero division 
    charfct = zeros(eltype(A),size(A))
    icsta = 1.0 - csta 
    iclta = 1.0 - clta 

    for ii in 2:N
        sq = A[ii] ^ 2 
        sta = csta * sq + icsta * sta 
        lta = clta * sq + iclta * lta 
        charfct[ii] = sta / lta 

        if ii < nlta 
            charfct[ii] = 0.0 
        end
    end 
    return charfct 
end

"""
  carl_sta_trig(A, nsta, nlta)

Computes the CarlSTAtrig characteristic function. 

eta = star - (ratio * ltar) - abs(sta - lta) - quiet

# Arguments
- `A::AbstractArray`: 1D array.
- `nsta:Int`: Length of short time average window in samples
- `nlta:Int`: Length of long time average window in samples
- `ratio:AbstractFloat`: as `ratio` gets smaller, carl_sta_trig gets more sensitive
- `Quiet:AbstractFloat`: as `quiet` gets smaller, carl_sta_trig gets more sensitive
"""
function carl_sta_trig(A::AbstractArray, nsta::Int, nlta::Int, ratio::Real, quiet::Real)
    N = size(A,1)
    T = eltype(A)
    sta = zeros(T,N)
    lta = deepcopy(sta)
    star = deepcopy(sta)
    ltar = deepcopy(sta)
    pad_sta = zeros(T,nsta)
    pad_lta = zeros(T,nlta)

    # compute the short time average (STA)
    for ii in 1:nsta
        sta .+= vcat(pad_sta,A[ii:N - nsta - 1 + ii])
    end
    sta ./= nsta

    # compute the long time average (LTA)
    for ii in 1:nlta
        lta .+= vcat(pad_lta, sta[ii:N - nlta - 1 + ii])
    end
    lta ./= nlta 
    lta = vcat(0,lta[1:end-1])
    
    # compute star, average of abs diff between trace and lta 
    for ii in 1:nsta
        star .+= vcat(pad_sta, abs.(A[ii:N - nsta - 1 + ii] .- lta[ii:N - nsta - 1 + ii]))
    end
    star ./= nsta 

    # compute ltar, long term average over star 
    for ii in 1:nlta 
        ltar .+= vcat(pad_lta, star[ii:N - nlta - 1 + ii])
    end 
    ltar ./= nlta 

    eta = star .- (ratio .* ltar) .- abs.(sta .- lta) .- quiet 
    eta[1:nlta] .= -1.0 
    return eta 
end

"""
  classic_sta_lta(A, nsta, nlta)

Computes classic STA/LTA. 

The length of the STA and LTA are given by `nsta` and `lsta` in samples, respectively.

# Arguments
- `A::AbstractArray`: 1D array.
- `nsta:Int`: Length of short time average window in samples
- `nlta:Int`: Length of long time average window in samples
"""
function classic_sta_lta(A::AbstractArray, nsta::Real, nlta::Real)
    T = eltype(A)
    # cumulative sum can be exploited to calculate a moving average 
    sta = cumsum(A .^ 2)
    lta = deepcopy(sta)

    # compute the STA and LTA 
    sta[nsta + 1:end] = sta[nsta + 1:end] .- sta[1:end-nsta]
    sta ./= nsta 
    lta[nlta + 1:end] = lta[nlta + 1:end] .- lta[1:end-nlta]
    lta ./= nlta 

    # pad zeros 
    sta[1:nlta] .= 0.0 

    # avoid division by zero by setting zero values to tiny float 
    tinyfloat = eps(T)
    ind = findall(lta .< tinyfloat)
    lta[ind] .= tinyfloat 
    return sta ./ lta 
end

"""
  delayed_sta_lta(A, nsta, nlta)

Delayed STA/LTA from Withers, 1998. 

# Arguments
- `A::AbstractArray`: 1D array.
- `nsta:Int`: Length of short time average window in samples
- `nlta:Int`: Length of long time average window in samples
"""
function delayed_sta_lta(A::AbstractArray, nsta::Int, nlta::Int)
    N = size(A,1)
    T = eltype(A)
    sta = zeros(T,N)
    lta = zeros(T,N)
    sta[1] = A[1] ^ 2 / nsta
    for ii in 2:N
        sta[ii] = (A[ii] ^ 2 + A[max(1,ii - nsta)] ^ 2) / nsta + sta[ii - 1]
        lta[ii] = (A[max(1,ii - nsta - 1)] ^ 2 + A[max(1, ii - nsta - nlta -1)] ^ 2) / nlta + lta[max(1,ii -1)]
    end
    sta[1:min(N,nlta + nsta + 50)] .= 0.0 
    lta[1:min(N,nlta + nsta + 50)] .= 1.0 # avoid division by zero 
    return sta ./ lta 
end

"""
  z_detect(A, nsta)

Z-detector from Withers, 1998. 

# Arguments
- `A::AbstractArray`: 1D array.
- `nsta:Int`: Length of z-detector time average window in samples
"""
function z_detect(A::AbstractArray, nsta::Int)
    N = size(A,1)
    T = eltype(A)

    # Z-detector given by Swindell and Snell, 1977 
    sta = zeros(T,N)
    pad_sta = zeros(T,nsta)

    for ii in 1:nsta 
        sta .+= vcat(pad_sta, A[ii:N - nsta - 1 + ii] .^ 2)
    end
    return (sta .- mean(sta)) ./ std(sta)
end

"""
  trigger_onset(charfct, thresh1, thresh2; max_len=1e99, max_len_delete=false)

Calculate trigger on and off times.

Given `thresh1` and `thresh2` calculate trigger on and off times from
characteristic function.
    
# Arguments
- `charfct::AbstractArray`: Characteristic function of e.g. STA/LTA trigger.
- `thresh1::Real`: Value above which trigger (of characteristic function)
    is activated (higher threshold).
- `thresh2::Real`: Value below which trigger (of characteristic function)
    is deactivated (lower threshold).
# Keywords
- `max_len::Int`: Maximum length of triggered event in samples. A new
    event will be triggered as soon as the signal reaches again above thresh1.
- `max_len_delete::Bool`: Do not write events longer than max_len into report file.
"""
function trigger_onset(charfct::AbstractArray, thresh1::Real, thresh2::Real;
    max_len::Int=10^10, max_len_delete::Bool=false,
)
    ind1 = findall(charfct .> thresh1)
    if length(ind1) == 0
        return []
    end
    ind2 = findall(charfct .> thresh2)

    on = [ind1[1]]
    off = [-1]

    # determine the indices where charfct falls below off-threshold
    ind2_ = Array{Bool}(undef,length(ind2))
    ind2_[1:end-1] .= diff(ind2) .> 1
    # last occurence is missed by diff, add it manually
    ind2_[end] = true
    append!(off,ind2[ind2_])
    append!(on,ind1[findall(diff(ind1) .> 1) .+ 1])
    # include last pick if trigger is on or drop it
    if max_len_delete
        # drop it
        append!(off,max_len)
        append!(on,on[end])
    else
        # include it
        append!(off,ind2[end])
    end

    pick = []
    while on[end] > off[1]
        while on[1] <= off[1]
            deleteat!(on,1)
        end
        while off[1] < on[1]
            deleteat!(off,1)
        end
        if off[1] - on[1] > max_len
            if max_len_delete
                deleteat!(on,1)
                continue
            end
            prepend!(off,on[1] + max_len)
        end
        push!(pick,[on[1],off[1]])
    end
    return permutedims(hcat(pick...))
end

function ar_pick(a, b, c, samp_rate, f1, f2, lta_p, sta_p, lta_s, sta_s, m_p, m_s,
    l_p, l_s, s_pick=true)
end

function coincidence_trigger()
end