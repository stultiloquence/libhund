# Future Work

- Add additional correctness tests
	- empty hyperedges
- Clean up build process
- Refactoring some overly complicated methods
	- especially the part where we combine permutations in run_multi_node ._.
- getting it to work on multiple nodes that are not on the devel partition
- Systematically run and evaluate and tweak using the tests
	- for different parameters / parameter variation to figure out a reasonable one
	- on cluster
- quality 2: relative increase in nnz (needs LU) -> compare with other libs
- performance test enhancement: profiling of different parts of the code
- Use external library for local reordering (currently returning identity permutations)
- Support additional matrix formats
- Unify the Loggger.gather() methods to all be void and then expose the result through getter function