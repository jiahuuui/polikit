#include <iostream>
#include <vector>

void poly_neighbor(int n_atom, std::vector<int> ref_n, std::vector<std::vector<int>> ref_list,
    std::vector<int> ptype, int x_type, int o_type, std::vector<int> &poly_n, std::vector<std::vector<int>> &poly_list) {
    for (int i = 0; i < n_atom; i++) {
        poly_n[i] = 0;
        poly_list[i].assign(20, 0);

        if (ptype[i] == o_type) {
            continue;
        }

        for (int n = 0; n < ref_n[i]; n++) {
            if (ptype[ref_list[i][n]] == o_type) {
                for (int m = 0; m < ref_n[ref_list[i][n]]; m++) {
                    if (ptype[ref_list[ref_list[i][n]][m]] == x_type && ref_list[ref_list[i][n]][m] != i) {
                        int k = 0;
                        while (true) {
                            if (poly_n[i] == 0) {
                                poly_n[i]++;
                                poly_list[i][0] = ref_list[ref_list[i][n]][m];
                                break;
                            }
                            if (ref_list[ref_list[i][n]][m] == poly_list[i][k]) {
                                break;
                            }
                            k++;
                            if (k >= poly_n[i]) {
                                poly_n[i]++;
                                poly_list[i][k] = ref_list[ref_list[i][n]][m];
                                if (poly_n[i] == 20) {
                                    std::cout << "Warning: Maximum neighbor length reached." << std::endl;
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
        // Sort the poly_list
        std::sort(poly_list[i].begin(), poly_list[i].begin() + poly_n[i]);
    }
    std::cout << "Polyhedra neighbor analysis -- done" << std::endl;
}

void neighbor_change(int n_atom, std::vector<int> old_n, std::vector<std::vector<int>> old_list,
    std::vector<int> new_n, std::vector<std::vector<int>> new_list, std::vector<int> &n_change) {
    for (int i = 0; i < n_atom; i++) {
        std::vector<int> unshared(old_n[i] + new_n[i]);
        int m = 0;
        int n = 0;
        int j = 0;
        while (m < old_n[i] && n < new_n[i]) {
            if (old_list[i][m] > new_list[i][n]) {
                unshared[j] = new_list[i][n];
                n++;
                j++;
            } else if (old_list[i][m] < new_list[i][n]) {
                unshared[j] = old_list[i][m];
                m++;
                j++;
            } else if (old_list[i][m] == new_list[i][n]) {
                m++;
                n++;
            }
        }
        while (m < old_n[i]) {
            unshared[j] = old_list[i][m];
            m++;
            j++;
        }
        while (n < new_n[i]) {
            unshared[j] = new_list[i][n];
            n++;
            j++;
        }
        j--;
        m = 0;
        n = 0;
        n_change[i] = j;
        n_change[i] = j;

        unshared.clear();
    }
}

void neighbor_change_ipr(int n_atom, std::vector<int> old_n, std::vector<std::vector<int>> old_list,
    std::vector<int> new_n, std::vector<std::vector<int>> new_list, std::vector<int> &n_change,
    std::vector<int64_t> &ipr_list, std::vector<double> &ipr) {
    for (int i = 0; i < n_atom; i++) {
        std::vector<int> unshared(old_n[i] + new_n[i]);
        int m = 0;
        int n = 0;
        int j = 0;
        int64_t nu = 0;
        int64_t de = 0;

        while (m < old_n[i] && n < new_n[i]) {
            if (old_list[i][m] > new_list[i][n]) {
                unshared[j] = new_list[i][n];
                n++;
                j++;
            } else if (old_list[i][m] < new_list[i][n]) {
                unshared[j] = old_list[i][m];
                m++;
                j++;
            } else if (old_list[i][m] == new_list[i][n]) {
                m++;
                n++;
                j++;
            }
        }

        while (m < old_n[i]) {
            unshared[j] = old_list[i][m];
            m++;
            j++;
        }

        while (n < new_n[i]) {
            unshared[j] = new_list[i][n];
            n++;
            j++;
        }

        j--;

        m = 0;
        n = 0;
        j = 0;
        n_change[i] = j;
        for (int n = 0; n < 200; n++) {
            nu += static_cast<int64_t>(ipr_list[i][n][1]) * static_cast<int64_t>(ipr_list[i][n][1]);
            de += static_cast<int64_t>(ipr_list[i][n][1]) * static_cast<int64_t>(ipr_list[i][n][1]);
        }

        if (de != 0) {
            ipr[i] = static_cast<double>(nu) / (static_cast<double>(de) * static_cast<double>(de));
        }
    }
}

void ln_analysis(int n_atom, std::vector<int> new_ln, std::vector<int> ref_ln, std::vector<int> nc, std::vector<std::vector<int>> &out) {
    for (int i = 0; i < n_atom; i++) {
        out[i][0] = 0;
        out[i][1] = 0;
        out[i][2] = 0;
        out[i][3] = 0;

        if (ref_ln[i] == 1 && new_ln[i] == 1) {
            out[i][0] += nc[i];
        } else if (ref_ln[i] == 1 && new_ln[i] == 2) {
            out[i][1] += nc[i];
        } else if (ref_ln[i] == 2 && new_ln[i] == 1) {
            out[i][2] += nc[i];
        } else if (ref_ln[i] == 2 && new_ln[i] == 2) {
            out[i][3] += nc[i];
        }
    }
}

void intersect(int n_atom, std::vector<int> poly_n, std::vector<std::vector<int>> poly_list,
    std::vector<int> ref_n, std::vector<std::vector<int>> ref_list, std::vector<int> &nout) {
    for (int i = 0; i < n_atom; i++) {
        nout[i] = 0;
        for (int j = 0; j < poly_n[i]; j++) {
            int k = 0;
            for (int m = 0; m < ref_n[i]; m++) {
                int n = 0;
                while (n < ref_n[poly_list[i][j]] && m < ref_n[i]) {
                    if (ref_list[i][m] > ref_list[poly_list[i][j]][n]) {
                        n++;
                    } else if (ref_list[i][m] < ref_list[poly_list[i][j]][n]) {
                        m++;
                    } else if (ref_list[i][m] == ref_list[poly_list[i][j]][n]) {
                        m++;
                        n++;
                        k++;
                    }
                }
            }
            if (nout[i] < k) {
                nout[i] = k;
            }
        }
    }
}

void classify(int n_atom, int out_length, std::vector<int> old_n, std::vector<int> new_n,
    std::vector<std::vector<int>> old_list, std::vector<std::vector<int>> new_list,
    int &n_increase, int &n_decrease, int &n_change, int n_start, int n_end, std::vector<int> &cn_out) {
    for (int j = 0; j < n_atom; j++) {
        if (old_n[j] < new_n[j]) {
            n_increase++;
        } else if (old_n[j] > new_n[j]) {
            n_decrease++;
        } else {
            for (int i = 0; i < old_n[j]; i++) {
                if (old_list[j][i] != new_list[j][i]) {
                    n_change++;
                    break;
                }
            }
        }
        for (int i = n_start; i <= n_end; i++) {
            if (new_n[j] == i) {
                cn_out[i - n_start]++;
                break;
            }
        }
    }
}
