# LatticeMC
Lattice Monte Carlo written in Julia

to perform Monte Carlo simulation on well defined surfaces (square, hexagonal lattice), then extract macroscopic rate information of reactions of interested.


author: Huijie Tian
email: hut216@lehigh.edu

## Installation
### Julia installation
	see https://julialang.org/downloads/

### Install LatticeMC
1. load LatticeMC
```
	git clone https://github.com/thj2009/LatticeMC.git
```
2. go to Julia REPL by
```
	julia
```

3. type "]" to go to package management, you should see something like this
```
	(@v1.6) pkg >
```

4. install LatticeMC and the dependencies
```
	add LatticeMC
```
5. test the installation with
```
	cd LatticeMC/test
	julia runtest.jl
```