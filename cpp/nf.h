// nf.h
#include <iostream>
#include <vector>

namespace polikit{
struct atom{
    int id;
    float x,y,z;
    int cn;
    std::vector<int> nlist;
};
class box{
public:
    int natom;
    float lx,ly,lz;
    float** xyz;
    float cutoff;
    void cvalue();
    void helpmessage();
    void neighbor_finder(double* xyz0, int atom_number, double lx, double ly, double lz, int pbc, double cutoff,
    int** neighbor, int* n_neighbor);

private:
};
}
