# PottsMC Julia Package

Monte Carlo simulation of Potts model.

## Installation

This code is written in Julia language. To add `PottsMC` package on Julia REPL run
```
Import("Pkg");Pkg.add("git@github.com:uguzzoni/PottsMC.git")
```

## Quick start

```julia
 
            using Random,PottsMC

            #random model
            A=10;L=5;
	        model = PottsMC.EpistasisModel(A,L); randn!(model.fields)

            #gradient descent
			s = rand(1:A,L); 
            PottsMC.descent!(s,model)

            #annealing with inverse temperature schedule
            betas=vcat(fill(0.1,100),fill(1.0,100),fill(10,100))
            PottsMC.anneal!(s,model,betas)

```

## License
[MIT license](LICENSE)
