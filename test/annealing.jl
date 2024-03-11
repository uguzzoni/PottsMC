using Test, Random, PottsMC

@testset "annealing" begin
	@testset "descent" begin
		for _=1:10
			A=rand(2:6); L=rand(2:6);
			
			model = PottsMC.FactorizedModel(A,L); randn!(model.fields)
			s = rand(1:A,L); PottsMC.descent!(s,model)
			for i=1:L, a=1:A
				st = copy(s); st[i] = a
				@test PottsMC.energy(model,st) ≥ PottsMC.energy(model,s)
			end

			model = PottsMC.EpistasisModel(A,L); randn!(model.fields)
			s = rand(1:A,L); PottsMC.descent!(s,model)
			for i=1:L, a=1:A
				st = copy(s); st[i] = a
				@test PottsMC.energy(model,st) ≥ PottsMC.energy(model,s)
			end
		end
	end
end