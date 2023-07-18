function generate(N)
	mapslices(rand(3,N); dims=1) do s
		 NamedTuple{(:cosθ,:ϕ,:Φ)}(s .* [2,2π,2π] .- [1,π,π])
	end |> vec |> DataFrame
end

arg(z) = atan(imag(z), real(z))
scale(t) = (1-cos(π*t))/2


function coordinates(px, basis) 
	p1, p2 = basis
	p = px .* cis(-arg(px[1])) #|> fold
	x = dot(p,p1) / sqrt(dot(p1,p1)*dot(p,p))
	y = dot(p,p2) / sqrt(dot(p2,p2)*dot(p,p))
	x, y
end

distance(a,b) = norm(a .- b)

function cluster(minima, ϵ=1e-5)
    selected = [minima[1]]
    for m in minima
        all(distance.(selected, Ref(m)) .> ϵ) && push!(selected, m)
    end
    return selected
end

function drop_conjugate(minima)
    return filter(minima) do m
        imag(m[2]) > 0
    end
end
