## 0. Plans and Logs

### 0.1 Logs

- *About r_cut and bin size* Because neighbor list construction and RDF analysis usually need very different bin capacity, I used an empirical formula to assign the bin capacity, which is `bin_cap = r_cut^3 * A`. Here `A` is an amplifier. We know the atomic density of alumina or galia is around 90/nm^3, so in average there can be `r_cut^3 * 90` atoms in a bin. For safety, we can have a constant factor of 20 for each bin, gives `r_cut^3 * 90 * 2` which can be used as bin_cap. But here the unit of `r_cut` is Angstrom not nm, by simplify and some approximation, we have `A = 2`, and the bin_cap is set as `bin_cap = r_cut^3 * 2`

- *About .xyz filetype* Data in each colume must be specified before analysis. Maybe a `.json` or `.toml` file is really needed to specify everything.

- *2024 Dec.* The RDF analysis can be performed now. The compute option is `g`, which means *g(r)* analysis. `.dump` and `.data` file can be read now.

3. find good way to dump atomic properties, in dynamic mode each static analysis dump the results in a text file.

4. decide how the comparison data is exported: for each dynamic analysis type write a specific comparison subroutine.

5. test script.

### 0.2 Plans

1. Pair-wise rdf. Structure factor S(k)

2. More flexible xyz/dump/data file import.

3. Atom type for pure silicon.

4. set capacity of neigh/poly_neigh as a constant.

### checklist:

poly
./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.32 -c p
bad
./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 2.32 -c b
rdf
./src/polikit -f ../test/ga2o3_test.xyz -p 1 -r 10 -c g
d_neigh_change
./src/polikit -d ../test/test_dir/ 3 -p 1 -r 2.3 -c p
## structure

```
parser.f90
    ^
    |
main.f90 ---> static.f90 -?-> nf.f90 -?-> ph.f90 -?-> bad.f90 -?-> rings.f90
    |            ^   |
[multi-frame?]   |   |-------------------------------|
    |            |                                   v
    --------> dynamic.f90 <-- compare.f90 <-- collect_data.f90
```

1. The `parser.f90` send back the compute option, dump option, file (directory) name to the `main.f90`.

2. The `data_input.f90` gets file name from `main.f90`, parse the file type and read the atomic coordination data. Data is stored in `data_input.f90`.

3. The static analysis modules (`neighbor_finder.f90, poly_analysis.f90, tct.f90, rings.f90, bad.f90, *etc*.`) will get atomic coordination data from `data_input.f90` and performed corresponding analysis. The results are stored in corresponding module, and at the end of the frame, they are destroyed.

4. The `dynamic.f90` is called when a directory instead of a file name is provided. It also includes the static analysis modules, and also has two more steps designed for a dynamic analysis, i.e. collect static analysis data and perform frame-wise comparison.
