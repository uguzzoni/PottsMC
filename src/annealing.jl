function mctrial!(s::AbstractVector{Int}, a::Integer,i::Integer,model::Model,β::Real)
	@assert β ≥ 0
	A = model.A; L = model.L
	@assert length(s) == L
	@assert eltype(s) <: Integer
	@assert 1 ≤ a ≤ A && 1 ≤ i ≤ L
	ΔE = Δenergy(model,s,a,i)
	if ΔE ≤ 0 || rand() < exp(-β*ΔE)
		s[i] = a
		return ΔE
	end
	return zero(ΔE)
end

function mc!(s::AbstractVector{Int},model::Model,β::Real)
	@assert β ≥ 0
	A = model.A; L = model.L
	@assert length(s) == L
	@assert eltype(s) <: Integer
	a = rand(1:A); i = rand(1:L);
	return mctrial!(s,a,i,model,β)
end

function descent!(s::AbstractVector{Int},model::Model)
	A = model.A; L = model.L
	localmin = false
	while !localmin
		localmin = true
		for a=1:A, i=1:L
			ΔE = Δenergy(model,s,a,i)
			if ΔE < 0
				s[i] = a
				localmin = false
				break
			end
		end
	end
end

function anneal!(seq::AbstractVector{Int},
                        model::PottsMC.Model,
                        schedule)
	A = model.A; L = model.L
    @assert length(seq) == L
    for i = 1:L @assert 1 ≤ seq[i] ≤ A end
    for β in schedule
        for i = 1 : length(seq)
            mc!(seq,model,β)
        end
    end
end