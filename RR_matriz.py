#! /usr/bin/env python

from numpy import linalg, zeros, dot, diag_indices_from, array, var, genfromtxt


def a_m(value, vector):
    mat = zeros((tr, tr))
    mat[diag_indices_from(mat)] = value
    mat = dot(dot(vector, mat), vector.transpose())
    return mat


def select_AV(eVal):
    """Selects the number of retained eigenvalues based on second
    derivative variance"""
    tr = len(eVal)
    eVal = eVal / sum(eVal)
    Bool = True
    grad = abs(eVal[:tr - 2] - 2 * eVal[1:tr - 1] + eVal[2:])
    Grad = open('./var_grad.txt', 'w')
    s1 = ''
    for i in range(2, tr - 5):
        aux = var(grad[i - 2:i + 3])
        s1 += str(tr - 2 - i) + ' ' + str(aux) + '\n'
        if(aux > (CutOff) and Bool):
            retained = tr - 2 - i
            Bool = False
    Grad.write(s1)
    Grad.flush()
    return retained


def NoiseControl(P, OutFile):
    """Extends Matrix P and prints it to OutFile"""
    eVal, eVec = linalg.eigh(P)
    tr = len(eVal)
    corte = select_AV(eVal)
    while  True:
        try:
            corte = int(raw_input("Type the number of eigenvalues to be\
 retained\nConsult the var_grad.txt file.\nSugestion: " + str(corte) + '\n'))
            break
        except ValueError:
            print "Invalid entry"
    eVal[:-corte] = eVal[corte]
    P = a_m(eVal, eVec)
    s1 = ''
    for i in range(tr):
        for j in range(tr):
            s1 += str(P[i, j]) + '\t'
        s1 += '\n'
    s1 += '\n'
    OutFile.write(s1)
    OutFile.flush()

CutOff = 1e-4

if __name__ == '__main__':
    #Number de Matrices#
    NP = 1
    #Input file#
    InputMatrix = open('./sam20.txt', 'r')
    Matrices = genfromtxt(InputMatrix)
    #Output File#
    Out = open('./matrix-non-noise.txt', 'w')

    tr = Matrices.shape[1]
    for k in range(NP):
        P = array(Matrices[k * tr:(k + 1) * tr])
        NoiseControl(P, Out)
