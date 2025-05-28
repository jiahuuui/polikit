# Polikit 0.4
[![GitHub](https://img.shields.io/badge/GitHub-V0.4-C71D23?logo=github&logoColor=white&labelColor=000)](https://github.com/jiahuuui/polikit/)
[![Bitbucket](https://img.shields.io/badge/Bitbucket-V0.4-0052CC?logo=bitbucket&logoColor=white&labelColor=000)](https://bitbucket.org/jiahuijiahui/polikit/src/master/)

## A polyhedral analysis toolkit

This package is originally developed for polyhedral analysis of amorphous structures. Now it has modules for other analysis methods, including bond angle analysis, RDF analysis, TCT analysis and ring statistics analysis. For now, file formats including xyz, lammps data file, lammps dump file can be read, but please carefully check the format when using the package. Analysis can be performed in either static or dynamic way, depends on whether the analysis only involves one file, or also comparison with other files.

Please contact **zjh239@foxmail.com** if you have bugs or issues to report.

#### 0.4  *(3 Mar. 2025)*

 - Non-affine displacement analysis.
 - Pair-wise cutoff values.
 - Dimension-wise periodic boundary condition.

#### 0.3  *(22 Jan. 2025)*

 - Ring statistics analysis.

#### 0.2
 - Bond angle distribution (BAD) analysis.
 - Radial distribution function.

#### 0.1
 - Polyhedral analysis.

## How to compile
At the root directory of the code:

1. `mkdir build && cd build`

2. `cmake ../.` or `cmake -DDEBUG=on ../.` for debug mode.

3. `make` or `make -j`

## What Polikit can do

### Analysis of a static configuration:

1. Polyhedral analysis

2. Bond angle distribution

3. Coordination distribution/change

4. Rings statistics distribution

5. Radial distribution function

6. Non-affine displacement

7. Cluster analysis

### Analysis of a series of configurations:

1. Topological constraint analysis

2. Neighbor change analysis

3. Polyhedral neighbor change analysis

## Usage
For static analysis(`-f`):

**`./polikit -f abc.xyz -p 1 -poly 2.3`**

`-p [int]` can be 1 or 0, decides whether periodic boundary condition will be applied.

`-[key] [parameter]` gives computing options and corresponding parameter. Now the availables are: `poly` - polyhedral analysis, `bad` - bond angle analysis, `rdf` - radial distribution function, `ring` - ring statistics analysis, `d2min` - non-affine displacement analysis, `cluster` - cluster analysis.

For dynamic analysis(`-d`):

**`./polikit -d dumpfiles -os 3 -p 1 -d2min 4.6`**

`-d` gives a directory name that contains .xyz files. 20 is the frame invertal for dynamic comparison. Other parameters work in the same way as in static analysis.

`-os [int]` indicates the frame interval, which is applicable when frame-wise comparison is performed.

## Examples

- Polyhedral analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -poly 2.3`

- Bond angle analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -bad 2.3`

- Radial distribution

`./polikit -f ../test/ga2o3_test.xyz -p 1 -rdf 10`

- Wendt-Abraham parameter calculation

`./polikit -f ../test/ga2o3_test.xyz -p 1 -wa 5`

- Honeycutt-Anderson parameters analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -ha 2.3`

- Ring statistics analysis

`./polikit -f ../test/ga2o3_test.xyz -p 1 -ring 2.3`

- Dynamic neighbor change analysis

`./polikit -d ../test/test_dir/ -os 1 -p 1 -nc 2.3`

- D2min analysis

`./polikit -d ../test/test_dir/ -os 2 -p 1 -d2min 4.6`

- D2min analysis and cluster analysis on high D2min atoms

`./polikit -d ../test/test_dir/ -os 2 -p 1 -d2min 4.6 -cluster 2.3`
