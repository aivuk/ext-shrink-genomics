#! /usr/bin/env python

from numpy import linalg, zeros, dot, diag_indices_from, array, var, genfromtxt

def a_m(value, vector):
    mat = zeros((tr, tr))
    mat[diag_indices_from(mat)] = value
    mat = dot(dot(vector, mat), vector.transpose())
    return mat

def last_big_eval_index(eig_vals, cutoff, var_matrix_filename):
    """Selects the number of retained eigenvalues based on second
    derivative variance"""
    len_ev = len(eig_vals)
    eig_vals_norm = eig_vals / sum(eig_vals)
    Bool = True
    second_diff = abs(eig_vals_norm[:len_ev - 2] - 2*eig_vals_norm[1:len_ev - 1] + eig_vals_norm[2:])
    variance_matrix_file = open(var_matrix_filename, 'w')
    s1 = ''

    for i in range(2, len_ev - 5):
        variance = var(second_diff[i - 2:i + 3])
        s1 += str(tr - 2 - i) + ' ' + str(variance) + '\n'
        if (variance > cutoff and Bool):
            last_index = tr - 2 - i
            Bool = False

    variance_matrix_file.write(s1)
    variance_matrix_file.flush()
    return last_index

def get_number(sug):
    ask_message = '''Type the number of eigenvalues to be
                retained\nConsult the var_grad.txt file.
                Sugestion:
                {}'''.format(sug)
    try:
        return int(raw_input(ask_message))
    except ValueError:
        print "Invalid entry"
        get_number(sug)

def noise_control(matrix, cutoff, output_filename, var_mat_filename):
    """Extends Matrix P and prints it to output_file"""
    eig_vals, eig_vecs = linalg.eigh(matrix)
    num_evals = len(eig_vals)
    last_ev_index = last_big_eval_index(eig_vals, cutoff, var_mat_filename)

    last_ev_index = get_number(last_ev_index)

    eig_vals[:-last_ev_index] = eig_vals[last_ev_index]
    p = a_m(eig_vals, eig_vecs)

    s1 = ''
    for i in range(num_evals):
        for j in range(num_evals):
            s1 += str(p[i, j]) + '\t'
        s1 += '\n'
    s1 += '\n'

    output_file = open(output_filename, 'w')
    output_file.write(s1)
    output_file.flush()

if __name__ == '__main__':
    #Number de Matrices#
    matrices_num = 1
    cutoff = 1e-4
    input_matrices = open('./sam20.txt', 'r')
    matrices = genfromtxt(input_matrices)
    new_matrix_filename = './matrix-non-noise.txt'
    var_matrix_filename = './var_grad.txt'

    tr = matrices.shape[1]
    for k in range(matrices_num):
        p = array(matrices[k * tr:(k + 1) * tr])
        noise_control(p, cutoff, new_matrix_filename,
                var_matrix_filename)
