# Polikit 0.2
## A polyhedral analysis toolkit

This package is originally developed for polyhedral analysis of amorphous structures. Now it has a module for analysis related to topological constraint theory also. It takes the dumped atomic configurations as input, and is able to perform static or dynamic analysis on the configurations.

#### 0.3 (working on ...)
 - Ring statistics analysis.

#### 0.2
 - Bond angle distribution (BAD) analysis.
 - Radial distribution function.

#### 0.1
 - Static TCT and polyhedral analysis.

## How to compile
At the root directory of the code:

1. `mkdir build && cd build`

2. `cmake ../.`

3. `make` or `make -j`

## What Polikit can do

### Analysis of a static configuration:

1. Polyhedral analysis

2. Bond angle distribution analysis

3. Coordination analysis

4. Rings statistics analysis

5. RDF

6. (Cavity analysis)

### Analysis of a series of configurations:

1. Topological constraint analysis

2. Neighbor change analysis

3. Polyhedral neighbor change analysis

## Usage
For static analysis(`-f`):

**`./polikit -f abc.xyz -p 1 -r 2.35 -c p`**

`-p [int]` can be 1 or 0, decides whether periodic boundary condition will be applied.

`-r [float]` indicates the cutoff distance.

`-c [string]` gives computing options, now the availables are: `p` - polyhedral analysis, `b` - bond angle analysis, `g` - radial distribution function.

For dynamic analysis(`-d`):

**`./polikit -d dumpfiles 20 -p 1 -r 2.35 -c pt`**

`-d` gives a directory name that contains .xyz files. 20 is the frame invertal for dynamic comparison. Other parameters work in the same way as in static analysis.

`-c` available computing options in dynamic analysis are: `t` - tct analysis.

## Test
  `./polikit -f ../../test_exmp/0.xyz -p 1 -r 2.35 -c pb`

  Perform polyhedral analysis on a single file.
<!--   **`./polikit -f ../../test_exmp/0.xyz -p 1 -r 2.35 -c t`** -->


  `./polikit -d ../../test_exmp/ -1 -p 1 -r 2.35 -c b`

  Perform BAD analysis for each file in a directory.

