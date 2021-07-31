# Persistence Calculator

A small but developing tool for fast computation of Persistence Diagram of Flag Filtrations. This work was inspired by [Ripser](https://github.com/Ripser/ripser) and uses some of the key ideas gathered there.

Here implemented basic algorithm of persistence diagram computation with some additional speed-ups inspired by other researches in computational topology:
 - Reducing the size of filtration by removing non-important simplices via the procedure of Edge Collapse (as it was [described](https://hal.inria.fr/hal-02873740) by J.-D. Boissonnat and S.Pritam). It's an optional preprocessing step 
 - Compute persistent *co*homology (as suggested by [Vin de Silva, Dmitriy Morozov, and Mikael Vejdemo-Johansson](https://doi.org/10.1088/0266-5611/27/12/124003))
 - Use the chain complex property that boundaries are cycles
    (employ the *clearing* optimization, aka *persistence with a twist*, as suggested by [Chao Chen and Michael Kerber](http://www.geometrie.tugraz.at/kerber/kerber_papers/ck-phcwat-11.pdf))
 - Parallelizng matrix reduction as suggested by [D. Morozov and A. Nigmetov.](https://www.mrzv.org/publications/lockfree-persistence/spaa/)


## Requirments:
This solution requires C++ template libraries [Boost](https://www.boost.org/) (version >= 1.71.0),  and [Intel Thread Building Blocks](http://www.github.com/oneapi-src/oneTBB)

## Building:
This project requires a C++17 compiler. Here is how this project can be obtained and built:
```sh
git clone https://github.com/ArGintum/PersistenceHomology.git
cd PersistenceHomology
make
```

## Running:
Example 
```sh
$ ./PersistenceCalculator input.txt [options]
```
Note that file

Available options:
 -  `--output`: name of file for output data (by default `output.txt` )
 -  `--format`: use the specified file format for the input. The following formats are supported:  
    - `matrix`: lower triangular distance matrix (default choice). Only space symbols as separators is acceptable at this version.
    - `cloud`: points cloud in a euclidean space. Only space symbols as separators is acceptable at this version.
 -  `--dim k`: Compute persistence homology up to dimension *k* (by defult 2)
 - `--threshold t`: compute Rips complexes up to diameter *t*
 - `--prep`: forced usage of Edge Collapse preprocessing (by default the necessity of preprocessing derived from input data and other parameters)
 - `--noprep`: forbid usage of Edge Collapse preprocessing (by default the necessity of preprocessing derived from input data and other parameters)
 - `--threads <n>`: number of threads for parallel computations (by default program uses all available CPU)

These options can be specified in any order.

Some examples:
```
$ ./PersistenceCalculator examples/H3N2_300.txt --output sample_output.txt --dim 4 --threads 3

...

$ ./PersistenceCalculator examples/H3N2_300.txt --dim 2 --threads 1 --noprep

...

$ ./PersistenceCalculator examples/H3N2_300.txt
```

## Plans for future (aka TODO list)

 - Improve performance of Edge Collapse procedure by changing model of graph representation
 - Add option of representative cycles for persistent homology
 - Add support of prime fields other than Z/*2*Z 
 - ...
 - Consider possible use of more "topological" way to simplify filtrations (maybe by exploring spectral sequences)


## Version

`1.0.3` July 2021.
