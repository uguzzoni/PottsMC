using Test, Random, LinearAlgebra, PottsMC

@testset "model" begin
    Random.seed!(1)
	@testset "energy" begin
		for _=1:10
			A=rand(2:6); L=rand(2:6); s = (rand(1:A,L)...,)
			a=rand(1:A); i=rand(1:L); st = ntuple(j -> j==i ? a : s[j],L)
	        
	        h = Dict((a,i)=>randn() for a=1:A for i=1:L)
	        Hh = -sum(h[(s[i],i)] for i=1:L)
	        fieldsh = [h[a,i] for i=1:L for a=1:A]
	        model = PottsMC.FactorizedModel(A,L,fieldsh)
	        @test PottsMC.energy(model,s) ≈ Hh
	        @test PottsMC.Δenergy(model,s,a,i) ≈ PottsMC.energy(model,st) - PottsMC.energy(model,s)
			harr = reshape(fieldsh, A, L) 
			@test harr == [h[a,i] for a=1:A, i=1:L]
			@test PottsMC.FactorizedModel(harr).fields == model.fields

	        J = Dict((a,b,i,j)=>randn() for a=1:A for b=1:A for i=1:L for j=i+1:L)
	        HJ = -sum(J[s[i],s[j],i,j] for j=2:L for i=1:j-1)
	        fieldsJ = [J[a,b,i,j] for j=2:L for i=1:j-1 for b=1:A for a=1:A]
	        model = PottsMC.EpistasisModel(A,L,[fieldsh;fieldsJ])
	        @test PottsMC.energy(model,s) ≈ HJ + Hh
	        @test PottsMC.Δenergy(model,s,a,i) ≈ PottsMC.energy(model,st) - PottsMC.energy(model,s)
			Jarr = [i == j ? 0.0 : i < j ? J[a,b,i,j]/2 : J[b,a,j,i]/2 for a=1:A, i=1:L, b=1:A, j=1:L]
			@test PottsMC.EpistasisModel(harr, Jarr).fields == model.fields
			
			# test energy is consistent with PhageMix
			S = rand(5:10)
			harr = randn(A,L)
			Jarr = randn(A,L,A,L)
			seqs = rand(1:A, L, S)
			s_oh = [seqs[i,s] == a for a=1:A, i=1:L, s=1:S]
			Epot = [PottsMC.energy(PottsMC.EpistasisModel(harr, Jarr), s) for s in eachcol(seqs)]
			Emix = -[dot(vec(s_oh[:,:,s]), reshape(Jarr, A*L, A*L), vec(s_oh[:,:,s])) for s=1:S]
			Emix -= reshape(s_oh, A * L, S)' * vec(harr)
			@test Emix ≈ Epot
		end
	end
end
