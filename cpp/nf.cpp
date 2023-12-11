// nf.cpp
// construct neighbor list from atom postions. Positions are read in import.cpp.
#include <iostream>
#include <vector>
#include <algorithm>
#include "nf.h"

using namespace std;
using namespace polikit;

void box::helpmessage(){
    cout<<100<<endl;
    cout<<"Example usage:     'x.exe -f abc.xyz -p 1 -r 2.32 -c pt -o ntpl'       (for single file analysis)"<<endl;
    cout<<"                   'x.exe -d xyzfiles 10 -p 1 -r 2.32 -c pt -o ntpl'   (for multi file analysis)"<<endl;
    cout<<"  -f [filename]    Load data which is a single file, incompatible with '-d' option."<<endl;
    cout<<"  -d [dir] [inv]   Load files in the directory, incompatible with '-f' option."<<endl;
    cout<<"  -r [float]       The cutoff value."<<endl;
    cout<<"  -p [1 .or. 0]    Apply periodic boundary condition or not."<<endl;
    cout<<"  -c [pt]          Compute option string. 'p' - polyhedral analysis; 't' - tct analysis."<<endl;
    cout<<"  -o [string]      Atomic output option string, won't dump atomic file if not set. 'n' - atomic coordination;"<<endl;
    cout<<"                   't' - tct results; 'p' - poly. neighbor; 'l' - linked state."<<endl;
// void function don not need return.
}

/*

void box::neighbor_finder(double* xyz0, int atom_number, double lx, double ly, double lz, int pbc, double cutoff,
    int** neighbor, int* n_neighbor) {

    // const int dp = 15; // Adjust as needed for selected_real_kind
    int o, p, q, xbin_max, ybin_max, zbin_max;
    vector<int> xbin(atom_number), ybin(atom_number), zbin(atom_number);
    bool isghost = false;
    double r, d;
    vector<std::vector<double>> realxyz(atom_number, std::vector<double>(3));
    vector<std::vector<double>> xyz(2 * atom_number, std::vector<double>(3));
    double x_min = *std::min_element(xyz0, xyz0 + atom_number);
    double y_min = *std::min_element(xyz0 + atom_number, xyz0 + 2 * atom_number);
    double z_min = *std::min_element(xyz0 + 2 * atom_number, xyz0 + 3 * atom_number);

    for (int i = 0; i < atom_number; i++) {
        realxyz[i][0] = xyz0[i] - x_min;
        realxyz[i][1] = xyz0[i + atom_number] - y_min;
        realxyz[i][2] = xyz0[i + 2 * atom_number] - z_min;
    }

    r = cutoff * cutoff;
    for (int i = 0; i < atom_number; i++) {
        xyz[i] = realxyz[i];
    }

    xbin_max = std::ceil(lx / cutoff) - 1;
    ybin_max = std::ceil(ly / cutoff) - 1;
    zbin_max = std::ceil(lz / cutoff) - 1;

    std::vector<std::vector<std::vector<bins>>>
        cellbins(xbin_max + 2, std::vector<std::vector<bins>>(ybin_max + 2, std::vector<bins>(zbin_max + 2)));

    for (int j = 0; j < xbin_max + 2; j++) {
        for (int i = 0; i < ybin_max + 2; i++) {
            for (int n = 0; n < zbin_max + 2; n++) {
                cellbins[j][i][n].n = 0;
                cellbins[j][i][n].n_list = std::vector<int>(500, 0);
            }
        }
    }

    int j = 1;
    std::vector<int> id(2 * atom_number), n, l, k;
    int xpbc = 0, ypbc = 0, zpbc = 0, xbin_t, ybin_t, zbin_t;

    for (int i = 0; i < atom_number; i++) {
        id[i] = i;
        xbin[i] = std::ceil(xyz[i][0] / cutoff);
        ybin[i] = std::ceil(xyz[i][1] / cutoff);
        zbin[i] = std::ceil(xyz[i][2] / cutoff);

        if (xbin[i] > xbin_max) {
            xbin[i] = xbin_max;
        }
        if (ybin[i] > ybin_max) {
            ybin[i] = ybin_max;
        }
        if (zbin[i] > zbin_max) {
            zbin[i] = zbin_max;
        }
        if (xbin[i] < 1) {
            xbin[i] = 1;
        }
        if (ybin[i] < 1) {
            ybin[i] = 1;
        }
        if (zbin[i] < 1) {
            zbin[i] = 1;
        }

        cellbins[xbin[i]][ybin[i]][zbin[i]].n++;
        cellbins[xbin[i]][ybin[i]][zbin[i]].n_list[cellbins[xbin[i]][ybin[i]][zbin[i]].n] = i;
        int n = 1;
        int l = 1;
        int k = 1;
        if (pbc == 1) {
            if (xbin[i] == 1) {
                xpbc = 1;
                isghost = true;
            }
            else if (xbin[i] == xbin_max) {
                xpbc = -1;
                n = -1;
                isghost = true;
            }
            if (ybin[i] == 1) {
                ypbc = 1;
                isghost = true;
            }
            else if (ybin[i] == ybin_max) {
                ypbc = -1;
                l = -1;
                isghost = true;
            }
            if (zbin[i] == 1) {
                zpbc = 1;
                isghost = true;
            }
            else if (zbin[i] == zbin_max) {
                zpbc = -1;
                k = -1;
                isghost = true;
            }
            if (isghost) {
                do {
                    for (int p = 0; p <= xpbc; p++) {
                        for (int q = 0; q <= ypbc; q++) {
                            for (int o = 0; o <= zpbc; o++) {
            if (p == 0 && q == 0 && o == 0) {
                continue;
            }
            id[atom_number + j] = i;
            xbin_t = xbin[i] + p * xbin_max;
            ybin_t = ybin[i] + q * ybin_max;
            zbin_t = zbin[i] + o * zbin_max;
            xyz[atom_number + j][0] = xyz[i][0] + lx * p;
            xyz[atom_number + j][1] = xyz[i][1] + ly * q;
            xyz[atom_number + j][2] = xyz[i][2] + lz * o;
            cellbins[xbin_t][ybin_t][zbin_t].n++;
            cellbins[xbin_t][ybin_t][zbin_t].n_list[cellbins[xbin[i]][ybin[i]][zbin[i]].n] = atom_number + j;
            j++;
                            }
                        }
                    }
                } while (false);
                xpbc = 0;
                ypbc = 0;
                zpbc = 0;
                isghost = false;
            }
        }
    }

    for (int i = 0; i < atom_number; i++) {
        n_neighbor[i] = 0;
        for (int j = 0; j < 10; j++) {
            neighbor[i][j] = 0;
        }
    }

    for (int i = 0; i < atom_number; i++) {
        n = 1;
        for (int p = -1; p <= 1; p++) {
            for (int q = -1; q <= 1; q++) {
                for (int o = -1; o <= 1; o++) {
                    if (cellbins[xbin[i] + p][ybin[i] + q][zbin[i] + o].n != 0) {
                        for (int k = 1; k <= cellbins[xbin[i] + p][ybin[i] + q][zbin[i] + o].n; k++) {
                            l = cellbins[xbin[i] + p][ybin[i] + q][zbin[i] + o].n_list[k];
                            if (i < id[l]) {
                                d = (xyz[i][0] - xyz[l][0]) * (xyz[i][0] - xyz[l][0]) +
                                    (xyz[i][1] - xyz[l][1]) * (xyz[i][1] - xyz[l][1]) +
                                    (xyz[i][2] - xyz[l][2]) * (xyz[i][2] - xyz[l][2]);
                                if (d < r) {
                                    n_neighbor[id[i]]++;
                                    n_neighbor[id[l]]++;
                                    neighbor[id[i]][n_neighbor[id[i]]] = id[l];
                                    neighbor[id[l]][n_neighbor[id[l]]] = id[i];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < atom_number; i++) {
        bubble_sort(n_neighbor[i], neighbor[i]);
    }

    // You can add print statements here for debugging if needed
}


void bubble_sort(int last, int array[]) {
    int temp;
    for (int i = last - 1; i >= 1; i--) {
        for (int j = 0; j < i; j++) {
            if (array[j + 1] < array[j]) {
                temp = array[j + 1];
                array[j + 1] = array[j];
                array[j] = temp;
            }
        }
    }
}*/
