# CoreNet
## Extended core characterization of a set of objects (metabolic reactions or genes) within a set of species

CoreNet performs probabilistic model based clusterings of objects (= metabolic reactions or genes) within a set of organisms. Two underlying model may be used: 1) a simple Bernoulli mixture model using only information about presence/absence of an object within a species ; 2) a hidden Markov random field model that combines two types of information: the former presence/absence information plus a dependency structure on the set of objects.
The dependency graph between reactions is the metabolic graph where two reactions are neighbors if they share a metabolite. The dependency graph between genes links any two genes whose positions are adjacent in the organism's genome.
