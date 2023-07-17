function NNL(m::Model, data)
	@unpack waveset = m
	I(τ) = intensity(m, τ)
	Nd = size(data,1)
	return -sum(log ∘ I, eachrow(data)) + Nd * sum(w->abs2(w.value), waveset)
end

function unfold(v::AbstractVector{<:Real})
    N = div(length(v)+1,2)
    return v[1:N] + 1im * [0, v[N+1:end]...]
end
function fold(v::AbstractVector{<:Number})
    [real(v)..., imag(v[2:end])...]
end

function go2min(m::Model, data, startvalues)
	objective(x) = NNL(update(m, x |> unfold), data)
	init = startvalues |> fold
	opt_result = optimize(objective, init, BFGS())
	return opt_result, opt_result.minimizer |> unfold 
end



struct Experiment{M,D,P,K}
    model::M
    data::D
    pv::SVector{K,SVector{P,C}} where C<:Number
end

nminima(exp::Experiment{M,D,P,K}) where {M,D,P,K} = K

function Experiment(m::Model, data, N::Int; initv::Vector{Vector{T}} where T)
    pickeditems = rand(1:size(data,1), N)
    _data = data[pickeditems,:]
    # 
    pv = [go2min(m, _data, init)[2] for init in initv]
    # 
    K = length(initv)
    P = length(initv[1])
    Experiment(m, _data, SVector{K}(SVector{P}.(pv)))
end

track(p0,pa,t) = (p0+pa)/2 + cis(π*t) * (p0-pa)/2


function NNL(exp::Experiment{M,D,P,K}, t) where {M,D,P,K}
    extended_pv = exp.pv[[K, 1:K..., 1]]
    itr_real = cubic_spline_interpolation(-1:K, extended_pv .|> real)
    itr_imag = cubic_spline_interpolation(-1:K, extended_pv .|> imag)
    # 
    f = itr_real(t) + 1im .* itr_imag(t)
    m = update(exp.model, f |> collect)
    return NNL(m, exp.data)
end
