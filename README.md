# Polikit 0.3
[![GitHub](https://img.shields.io/badge/GitHub-V0.3-C71D23?logo=github&logoColor=white&labelColor=000)](https://github.com/jiahuuui/polikit/)
[![Bitbucket](https://img.shields.io/badge/Bitbucket-V0.3-0052CC?logo=bitbucket&logoColor=white&labelColor=000)](https://bitbucket.org/jiahuijiahui/polikit/src/master/)
<!-- [![GitLab](https://img.shields.io/badge/GitLab-Repository-FFD700?logo=gitlab&logoColor=white&labelColor=DC143C)](https://gitlab.com/jhcheung/polikit) [![Gitee](https://img.shields.io/badge/Gitee-Repository-FFD700?logo=gitee&logoColor=white&labelColor=DC143C)](https://gitlab.com/jiahuiiii)-->

## A polyhedral analysis toolkit

This package is originally developed for polyhedral analysis of amorphous structures. Now it has a module for analysis related to topological constraint theory also. It takes the dumped atomic configurations as input, and is able to perform static or dynamic analysis on the configurations.

#### 0.3
  *(Released on 22 Jan. 2025)*
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


### Analysis of a series of configurations:

1. Topological constraint analysis

2. Neighbor change analysis

3. Polyhedral neighbor change analysis

## Usage
For static analysis(`-f`):

**`./polikit -f abc.xyz -p 1 -r 2.35 -c pbgr`**

`-p [int]` can be 1 or 0, decides whether periodic boundary condition will be applied.

`-r [float]` indicates the cutoff distance.

`-c [string]` gives computing options, now the availables are: `p` - polyhedral analysis, `b` - bond angle analysis, `g` - radial distribution function, `r` - ring statistics analysis.

For dynamic analysis(`-d`):

**`./polikit -d dumpfiles 20 -p 1 -r 2.35 -c pt`**

`-d` gives a directory name that contains .xyz files. 20 is the frame invertal for dynamic comparison. Other parameters work in the same way as in static analysis.

`-c` available computing options in dynamic analysis are: `t` - tct analysis.

## Examples

- Polyhedral analysis

`./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.32 -c p`

- Bond angle analysis

`./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.32 -c b`

- Radial distribution

`./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 10 -c g`

- Ring statistics analysis

`./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.3 -c r`

- Dynamic neighbor change analysis

`./src/polikit -d ../test/test_dir/ 3 -p 1 -r 2.3 -c p`
