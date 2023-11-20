#include <iostream>
#include <cmath>
#include <vector>

// Function to calculate standard deviation from a 1D array
template<typename T>
T stddev(const std::vector<T>& array) {
    T bar = 0;
    T sm = 0;
    for (const T& value : array) {
        bar += value;
    }
    int q = array.size();
    bar /= q;
    for (const T& value : array) {
        sm += std::pow(value - bar, 2);
    }
    return std::sqrt(sm / q);
}

// Function to convert stable angle number to ncbb
template<typename T>
T angle2ncbb(T angle_number) {
    if (angle_number < 1.0) {
        return 0.0;
    } else if (angle_number < 3.0) {
        return 1.0;
    } else if (angle_number < 6.0) {
        return 3.0;
    } else if (angle_number < 10.0) {
        return 5.0;
    } else if (angle_number < 15.0) {
        return 7.0;
    } else if (angle_number < 21.0) {
        return 9.0;
    } else if (angle_number < 28.0) {
        return 11.0;
    } else {
        return 13.0;
    }
}

void tct_calculate(int n_atom, const std::vector<std::vector<std::vector<double>>&> vectors,
    std::vector<std::vector<std::vector<double>>&> bond_length,
    std::vector<std::vector<std::vector<double>>&> bond_angle, const std::vector<int>& n_neighbor,
    const std::vector<std::vector<int>>& neighbor, std::vector<std::vector<int>>& length_header,
    std::vector<std::vector<std::vector<int>>&> angle_header, std::vector<std::vector<std::vector<int>>&> counter,
    std::vector<std::vector<double>>& constrain) {

    const double pi = 3.141592653589793;
    const double r2d = 180.0 / pi;

    for (int i = 0; i < n_atom; i++) {
        for (int l = 19; l >= 2; l--) {
            for (int j = 0; j < 10; j++) {
                bond_length[i][j][l] = bond_length[i][j][l - 1];
            }
        }

        for (int l = 19; l >= 2; l--) {
            for (int j = 0; j < 30; j++) {
                bond_angle[i][j][l] = bond_angle[i][j][l - 1];
            }
        }

        for (int j = 0; j < 10; j++) {
            if (counter[i][j][1] == 0 && length_header[i][j] > 0) {
                for (int l = j; l < 9; l++) {
                    for (int k = 0; k < 2; k++) {
                        if (angle_header[i][l][k] == length_header[i][j]) {
                            for (int m = 0; m < 30; m++) {
                                angle_header[i][l][k] = angle_header[i][l + 1][k];
                                bond_angle[i][l][k] = bond_angle[i][l + 1][k];
                            }
                            if (l < 29) {
                                l--;
                            }
                            break;
                        }
                    }
                }
                length_header[i][j] = 0;
                for (int l = j; l < 9; l++) {
                    for (int k = 0; k < 2; k++) {
                        angle_header[i][l][k] = angle_header[i][l + 1][k];
                        bond_angle[i][l][k] = bond_angle[i][l + 1][k];
                    }
                }
                angle_header[i][29][0] = 0;
                angle_header[i][29][1] = 0;
                bond_angle[i][29][0] = 0;
                bond_angle[i][29][1] = 0;
                j--;
            }
        }

        for (int j = 0; j < 10; j++) {
            if (counter[i][j][1] < 20) {
                counter[i][j][1] = 0;
            }
            if (counter[i][j][2] > 0) {
                counter[i][j][2]--;
            }
        }

        for (int j = 0; j < 10; j++) {
            int length = neighbor[i][j];
            double length_val = 0.0;
            for (int m = 0; m < 3; m++) {
                length_val += vectors[i][j][m] * vectors[i][j][m];
            }
            length_val = std::sqrt(length_val);
            for (int l = 0; l < 10; l++) {
                if (neighbor[i][j] == length_header[i][l]) {
                    bond_length[i][j][0] = length_val;
                    counter[i][j][0]++;
                    counter[i][j][1] = 20;
                    break;
                }
                else if (neighbor[i][j] < length_header[i][l]) {
                    for (int k = 9; k >= l; k--) {
                        length_header[i][k + 1] = length_header[i][k];
                        for (int m = 0; m < 20; m++) {
                            bond_length[i][j][m + 1] = bond_length[i][j][m];
                        }
                    }
                    length_header[i][l] = neighbor[i][j];
                    bond_length[i][j][0] = length_val;
                    counter[i][j][0] = 1;
                    counter[i][j][1] = 20;
                    break;
                }
                else if (length_header[i][l] == 0) {
                    length_header[i][l] = neighbor[i][j];
                    bond_length[i][j][0] = length_val;
                    counter[i][j][0]++;
                    counter[i][j][1] = 20;
                    break;
                }
            }
        }

        for (int j = 0; j < n_neighbor[i]; j++) {
            double length_a = 0.0;
            for (int m = 0; m < 3; m++) {
                length_a += vectors[i][j][m] * vectors[i][j][m];
            }
            length_a = std::sqrt(length_a);

            for (int l = j + 1; l < n_neighbor[i]; l++) {
                double angle = r2d * std::acos(
                    (vectors[i][j][0] * vectors[i][l][0] + vectors[i][j][1] * vectors[i][l][1] +
                    vectors[i][j][2] * vectors[i][l][2]) /
                    (length_a * norm2(vectors[i][l])));

                for (int k = 0; k < 30; k++) {
                    if (length_header[i][k] == length_a && angle_header[i][k][0] == neighbor[i][j] && angle_header[i][k][1] == neighbor[i][l]) {
                        bond_angle[i][k][0] = angle;
                        break;
                    }
                    else if (length_a < length_header[i][k]) {
                        for (int m = 29; m >= k; m--) {
                            angle_header[i][m + 1][0] = angle_header[i][m][0];
                            angle_header[i][m + 1][1] = angle_header[i][m][1];
                            for (int n = 0; n < 20; n++) {
                                bond_angle[i][m + 1][n] = bond_angle[i][m][n];
                            }
                        }
                        angle_header[i][k][0] = length_a;
                        angle_header[i][k][1] = neighbor[i][j];
                        bond_angle[i][k][0] = angle;
                        break;
                    }
                    else if (length_header[i][k] == 0) {
                        angle_header[i][k][0] = length_a;
                        angle_header[i][k][1] = neighbor[i][j];
                        bond_angle[i][k][0] = angle;
                        break;
                    }
                }
            }
        }

        int tm = 0;
        int tm_list[10] = {0};

        for (int j = 0; j < 10; j++) {
            if (counter[i][j][0] >= 20 && counter[i][j][1] == 20) {
                tm++;
                tm_list[tm - 1] = length_header[i][j];
                double tct = stddev(bond_length[i][j]);
                if (tct < 1.0) {
                    constrain[i][0] += 1.0;
                }
                for (int l = j + 1; l < 10; l++) {
                    if (counter[i][l][0] >= 20 && counter[i][l][1] == 20) {
                        for (int k = 0; k < 30; k++) {
                            if (angle_header[i][k][0] == length_header[i][j] && angle_header[i][k][1] == length_header[i][l]) {
                                double tct = stddev(bond_angle[i][k]);
                                if (tct < 15.0) {
                                    constrain[i][1] += 1.0;
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }

        constrain[i][0] = 0.5 * constrain[i][0];

        if (constrain[i][1] >= 1.0) {
            constrain[i][1] = angle2ncbb(constrain[i][1]);
        }
    }
}

// int main() {
//     // Example usage of tct_calculate function with appropriate data structures
//     int n_atom = 100; // Replace with your actual number of atoms
//     std::vector<std::vector<std::vector<double>>> vectors(n_atom, std::vector<std::vector<double>>(10, std::vector<double>(3, 0.0)));
//     std::vector<std::vector<std::vector<double>>> bond_length(n_atom, std::vector<std::vector<double>>(10, std::vector<double>(20, 0.0)));
//     std::vector<std::vector<std::vector<double>>> bond_angle(n_atom, std::vector<std::vector<double>>(30, std::vector<double>(20, 0.0)));
//     std::vector<int> n_neighbor(n_atom, 10); // Replace with your actual neighbor data
//     std::vector<std::vector<int>> neighbor(n_atom, std::vector<int>(10, 0)); // Replace with your actual neighbor data
//     std::vector<std::vector<int>> length_header(n_atom, std::vector<int>(10, 0));
//     std::vector<std::vector<std::vector<int>>> angle_header(n_atom, std::vector<std::vector<int>>(30, std::vector<int>(2, 0)));
//     std::vector<std::vector<std::vector<int>>> counter(n_atom, std::vector<std::vector<int>>(10, std::vector<int>(2, 0)));
//     std::vector<std::vector<double>> constrain(n_atom, std::vector<double>(2, 0.0));
//
//     tct_calculate(n_atom, vectors, bond_length, bond_angle, n_neighbor, neighbor, length_header, angle_header, counter, constrain);
//
//     // Output the results in the 'constrain' vector
//     for (int i = 0; i < n_atom; i++) {
//         std::cout << "Atom " << i << ": Constraint 1 = " << constrain[i][0] << ", Constraint 2 = " << constrain[i][1] << std::endl;
//     }
//
//     return 0;
// }

void poly_volume(int n_atom, std::vector<int>& n_neigh, std::vector<std::vector<std::vector<double>>>& delta, int ctype, int n_ctype, std::vector<double>& hull_v, std::vector<double>& std_volume) {
    int dp = 15;
    std::vector<std::vector<double>> edge(n_atom, std::vector<double>(20, 0.0));

    double vfactor = 6.0 * std::sqrt(2.0);
    int c = 0;

    for (int n = 0; n < n_atom; ++n) {
        if (ptype[n] == ctype && n_neigh[n] == 4) {
            int k = 0;
            for (int i = 0; i < n_neigh[n] - 1; ++i) {
                for (int j = i + 1; j < n_neigh[n]; ++j) {
                    double norm = 0.0;
                    for (int d = 0; d < 3; ++d) {
                        double diff = delta[n][j][d] - delta[n][i][d];
                        norm += diff * diff;
                    }
                    edge[n][k] = std::sqrt(norm);
                    k++;
                }
            }

            double mean = 0.0;
            for (int i = 0; i < k; ++i) {
                mean += edge[n][i];
            }
            mean /= k;

            std_volume[c] = hull_v[c] / std::pow(mean, 3) * vfactor;
            c++;
        }
    }
}

// int main() {
//     // Example usage
//     int n_atom = 10;
//     std::vector<int> n_neigh(n_atom, 4);
//     std::vector<std::vector<std::vector<double>>> delta(n_atom, std::vector<std::vector<double>>(10, std::vector<double>(3, 0.0)));
//     int ctype = 1;
//     int n_ctype = 2;
//     std::vector<double> hull_v(n_ctype, 0.0);
//     std::vector<double> std_volume(n_ctype, 0.0);
//
//     poly_volume(n_atom, n_neigh, delta, ctype, n_ctype, hull_v, std_volume);
//
//     return 0;
// }
