#! /usr/bin/env python
"""Usage: ./rr_matrix.py [options]

Options:
    -h --help                     show this      
    -c CUTOFF, --cutoff=CUTOFF    cuttoff value [default: 1e-4]
    -i input_file                 input matrices filename [default: sam20.txt]
    -o output_file                output file [default: matrix-non-noise.txt]
    -v variance_file              variance matrix file [default: var_grad.txt]
"""

from numpy import linalg, zeros, dot, diag_indices_from, array, var, genfromtxt
from docopt import docopt
    
def a_m(value, vector):
    mat = zeros((len(value), len(value)))
    mat[diag_indices_from(mat)] = value
    mat = dot(dot(vector, mat), vector.transpose())
    return mat

def last_big_eval_index(eig_vals, cutoff, var_matrix_filename):
    """Selects the number of retained eigenvalues based on second
    derivative variance"""
    len_ev = len(eig_vals)
    eig_vals_norm = eig_vals / sum(eig_vals)
    second_diff = abs(eig_vals_norm[:len_ev - 2] - 2*eig_vals_norm[1:len_ev - 1] + eig_vals_norm[2:])
    variance_matrix_file = open(var_matrix_filename, 'w')
    vars_txt = ''
    last_index = None

    for i in range(2, len_ev - 5):
        variance = var(second_diff[i - 2:i + 3])
        vars_txt += "{} {}\n".format(len_ev - 2 - i, variance)
        if not last_index and variance > cutoff:
            last_index = len_ev - 2 - i

    variance_matrix_file.write(vars_txt)
    variance_matrix_file.flush()
    return last_index

def ask_number(sug):
    ask_message = """
    Type the number of eigenvalues to be
    retained. 
    Consult the var_grad.txt file.
    Sugestion: {} > """.format(sug)
    try:
        return int(raw_input(ask_message))
    except ValueError:
        print "\nInvalid entry"
        ask_number(sug)

def noise_control(matrix, cutoff, output_filename, var_mat_filename):
    """Extends Matrix P and prints it to output_file"""
    eig_vals, eig_vecs = linalg.eigh(matrix)
    num_evals = len(eig_vals)
    last_ev_index = last_big_eval_index(eig_vals, cutoff, var_mat_filename)
    last_ev_index = ask_number(last_ev_index)
    eig_vals[:-last_ev_index] = eig_vals[last_ev_index]
    p = a_m(eig_vals, eig_vecs)
    
    new_matrix = ''
    for i in range(num_evals):
        for j in range(num_evals):
            new_matrix += "{}\t".format(p[i, j])
        new_matrix += '\n'
    new_matrix += '\n'

    output_file = open(output_filename, 'w')
    output_file.write(new_matrix)
    output_file.flush()

def main(options):
    # Number de Matrices
    matrices_num = 1
    cutoff = float(options['--cutoff'])
    input_matrices = open(options['-i'], 'r')
    matrices = genfromtxt(input_matrices)
    new_matrix_filename = options['-o'] 
    var_matrix_filename = options['-v']  

    tr = matrices.shape[1]
    for k in range(matrices_num):
        p = array(matrices[k * tr:(k + 1) * tr])
        noise_control(p, cutoff, new_matrix_filename,
                var_matrix_filename) 

if __name__ == '__main__':
    options = docopt(__doc__)
    main(options)
