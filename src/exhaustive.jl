function generate_sequences(A::Integer,L::Integer)
    @assert L ≥ 0 && A ≥ 1
    Iterators.flatten((Iterators.product((1:A for i in 1:L)...),))
end
collect_sequences(A::Integer,L::Integer) = collect(generate_sequences(A,L))

function exhaustive_groundstateFactorized(A::Integer,L::Integer,fields)
	@assert length(fields) == A*L
	s0 = ntuple(i->0,L)
	for a=1:A, i=1:L
		fields[a + (i-1)*A]
	end
end