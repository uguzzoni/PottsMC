abstract type Model{U<:Real} end
# A type parameter U allows us to use ForwardDiff (which sends a DualNumber<:Real)

struct EpistasisModel{U} <: Model{U}
    A::Int; L::Int; fields::Vector{U};
    function EpistasisModel{U}(A::Integer,L::Integer,fields::AbstractVector{U}=zeros(U,fieldslen(EpistasisModel,A,L))) where {U}
        @assert A ≥ 0 && L ≥ 0
        @assert length(fields) == fieldslen(EpistasisModel,A,L)
        new(A,L,fields)
    end
end
EpistasisModel(A::Integer,L::Integer,fields::AbstractVector{U}) where {U<:Real} = EpistasisModel{U}(A,L,fields)
EpistasisModel(A::Integer,L::Integer) = EpistasisModel{Float64}(A,L)

function EpistasisModel(h::AbstractArray{<:Real,2}, J::AbstractArray{<:Real,4})
	A, L = size(h)
	@assert size(J) == (A,L,A,L)
	h_ = vec(h) + [J[a,i,a,i] for i = 1:L for a = 1:A]
	J_ = [J[a,i,b,j] + J[b,j,a,i] for j=2:L for i=1:j-1 for b=1:A for a=1:A]
	return EpistasisModel(A, L, [h_; J_])
end

struct FactorizedModel{U} <: Model{U}
    A::Int; L::Int; fields::Vector{U};
    function FactorizedModel{U}(A::Integer,L::Integer,fields=zeros(U,fieldslen(FactorizedModel,A,L))) where {U}
        @assert A ≥ 0 && L ≥ 0
        @assert length(fields) == fieldslen(FactorizedModel,A,L)
        new(A,L,fields)
    end
end
FactorizedModel(A::Integer,L::Integer,fields::AbstractVector{U}) where {U<:Real} = FactorizedModel{U}(A,L,fields)
FactorizedModel(A::Integer,L::Integer) = FactorizedModel{Float64}(A,L)

function FactorizedModel(h::AbstractMatrix{<:Real})
	A, L = size(h)
	return FactorizedModel(A, L, vec(h))
end

fieldslen(::Type{EpistasisModel},A::Integer,L::Integer) = A*(L + A * (((L-1) * L) >> 1))
fieldslen(::Type{FactorizedModel},A::Integer,L::Integer) = A*L
fieldslen(::Type{EpistasisModel{U}},A::Integer,L::Integer) where {U} = fieldslen(EpistasisModel,A,L)
fieldslen(::Type{FactorizedModel{U}},A::Integer,L::Integer) where {U} = fieldslen(FactorizedModel,A,L)

function energy(model::FactorizedModel{U},sequence) where {U}
	@assert eltype(sequence) <: Integer
	@assert length(sequence) == model.L
	A = model.A; L = model.L; H = zero(U)
	for i=1:L
		@assert 1 ≤ sequence[i] ≤ A
		#= I think that the CPU is able to do (i-1)*A as fast as a single 
		integer addition, so no need to have an independent 'offset' variable
		here. However we might want to benchmark a few variations here  =#
		H -= model.fields[sequence[i] + (i-1)*A] 
	end
	return H
end
function energy(model::EpistasisModel{U},sequence) where {U}
	A = model.A; L = model.L; A2 = A^2
	@assert eltype(sequence) <: Integer 
	@assert length(sequence) == L
	H = zero(U); offsetJ = A*L
	for j=1:L
		@assert 1 ≤ sequence[j] ≤ A
		H -= model.fields[sequence[j] + (j-1)*A]
		local_offset = (sequence[j]-1)*A
		for i = 1:j-1
			H -= model.fields[sequence[i] + local_offset + offsetJ]
			offsetJ += A2
		end
	end
	return H
end

"ΔE of flipping site i to color a"
function Δenergy(model::FactorizedModel{U},sequence,a::Integer,i::Integer) where {U}
	A = model.A; L = model.L; A2 = A^2
	@assert eltype(sequence) <: Integer 
	@assert length(sequence) == L
	@assert 1 ≤ a ≤ A
	@assert 1 ≤ i ≤ L
	@assert 1 ≤ sequence[i] ≤ A
	offset = (i-1)*A
	ΔE = model.fields[sequence[i] + offset] - model.fields[a + offset]
	return ΔE

end
function Δenergy(model::EpistasisModel{U},sequence,a::Integer,i::Integer) where {U}
	A = model.A; L = model.L; A2 = A^2; AL = A*L
	@assert eltype(sequence) <: Integer 
	@assert length(sequence) == L
	@assert 1 ≤ a ≤ A
	@assert 1 ≤ i ≤ L
	@assert 1 ≤ sequence[i] ≤ A
	
	offset = (i-1)*A
	ΔE = model.fields[sequence[i] + offset] - model.fields[a + offset]

	offset = AL + binomial(i-1,2)*A2
	for j=1:i-1
		@assert 1 ≤ sequence[j] ≤ A
		local_offset = offset + sequence[j] + (j-1)*A2
		ΔE += model.fields[local_offset + (sequence[i]-1)*A] - model.fields[local_offset + (a-1)*A]
		#ΔE += model.fields[AL + sequence[j] + (sequence[i]-1)*A + (j-1)*A2 + binomial(i-1,2)*A2] - model.fields[AL + sequence[j] + (a-1)*A + (j-1)*A2 + binomial(i-1,2)*A2]
	end
	offset += (i-1)*A2
	for j=i+1:L
		@assert 1 ≤ sequence[j] ≤ A
		offset += (j-2)*A2
		local_offset = offset + (sequence[j]-1)*A
		ΔE += model.fields[local_offset + sequence[i]] - model.fields[local_offset + a]
		#ΔE += model.fields[AL + sequence[i] + (sequence[j]-1)*A + (i-1)*A2 + binomial(j-1,2)*A2] - model.fields[AL + a + (sequence[j]-1)*A + (i-1)*A2 + binomial(j-1,2)*A2]
	end

	return ΔE
end