'''Implementation of
An Explicit Algorithm for Computing the Picard Number of Certain Algebraic Surfaces
Author(s): Tetsuji Shioda
https://www.jstor.org/stable/pdf/2374678.pdf?refreqid=excelsior%3Aa5ea099e65db9eae2ba9c2d655762b4c'''

import numpy as np
import itertools

from util.compositions import compositions
from log import get_logger
import logging

log = get_logger(__name__, with_logfile=True, level=logging.INFO)


def matrix_cofactor(matrix):
    '''
    Computes the cofactor matrix. Given a matrix A, its cofactor matrix B is the
    unique matrix satisfying AB=det(A)*Id.
    :param matrix: An array-like quadratic matrix
    :return: Its cofactor matrix
    '''
    return np.linalg.inv(matrix) * np.linalg.det(matrix)

def get_delta(matrix):
    '''
    Compute the integer delta defined on page 417 of the article
    "An Explicit Algorithm ..."
    '''
    cof = np.rint(matrix_cofactor(matrix)).astype('int16')
    delta = np.gcd.reduce(cof.flatten())
    return delta

def get_d(matrix):
    '''
    Compute the integer d defined on page 417 of the article
    "An Explicit Algorithm ..."
    '''
    d = abs(np.linalg.det(matrix))/get_delta(matrix)
    return int(np.round(d))

def get_B_matrix(matrix):
    '''
    Compute the module B^2_d defined on page 418 of the article
    "An Explicit Algorithm ..."
    '''
    d = get_d(matrix)
    B = np.rint(d*np.linalg.inv(matrix)).astype('int16')
    return B

def get_iterator_of_A():
    '''
    Gives a generator for all 4x4-matrices with the property that the
    entries in each line sum up to 4. Matrices need not define Delsarte
    surfaces, because they may be singular, for example.
    :return: A generator for all 4x4-matrices with the property that the
    entries in each line sum up to 4
    '''
    prod_generator = itertools.product(compositions(4, 4),compositions(4, 4),
                                       compositions(4, 4),compositions(4, 4))
    return prod_generator

def check_valid_A(matrix):
    '''
    Checks if a 4x4 matrix A in which all rows sum to 4 defines a smooth Delsarte
    K3 surface. It is checked that
    (1) det(A)=0, (2) each column contains a zero, (3) each column contains a
    3 or 4, which is equivalent to smoothness. The condition that each row sums
    up to 4 is not checked.
    :param matrix: A 4x4 positive integer matrix in whichall rows sum to 4
    :return: True if the matrix defines a smooth Delsarte K3 surface, False otherwise
    '''
    if np.rint(np.linalg.det(matrix)) == 0:
        return False
    if not np.array_equal(matrix.min(axis=0), [0,0,0,0]):
        return False

    # Smoothness check: check if each column contains a 3 or a 4.
    # For example, if first column contains no 3 and no 4, then
    # [1:0:0:0] is singular point
    if np.all(np.greater_equal(matrix.max(axis=0), [3,3,3,3])):
        return True
    else:
        return False

def get_M_d_module(matrix):
    '''
    Compute the module M_d defined on page 417 of the article
    "An Explicit Algorithm ..."
    '''
    M = set()
    d = get_d(matrix)
    for tuple in itertools.product(range(d), repeat=3):
        M.add((tuple[0], tuple[1], tuple[2], (-tuple[0]-tuple[1]-tuple[2])%d))
    return M

def get_L_A_module(matrix):
    '''
    Compute the module L_A defined on page 417 of the article
    "An Explicit Algorithm ..."
    '''
    L = set()
    B = get_B_matrix(matrix)
    d = get_d(matrix)
    for tup in get_M_d_module(matrix):
        new_el = list(np.mod(np.matmul(np.array(tup), B), d))
        L.add(tuple(new_el))
    return L

def get_A_squared_d(matrix):
    '''
    Compute the module A^2_d defined on page 417 of the article
    "An Explicit Algorithm ..."
    '''
    A = set()
    log.debug('d=%s, delta=%s' % (get_d(matrix), get_delta(matrix)))
    M = get_M_d_module(matrix)
    log.debug('len(M)=%s' % len(M))
    for tup in M:
        if tup[0] != 0 \
                and tup[1] != 0 \
                and tup[2] != 0 \
                and tup[3] != 0:
            A.add(tup)
    log.debug('len(A)=%s' % len(A))
    return A

def get_B_squared_d(matrix):
    '''
    Compute the module B^2_d defined on page 417 of the article
    "An Explicit Algorithm ..."
    '''
    B = set()
    d = get_d(matrix)
    A = get_A_squared_d(matrix)
    for (k,tup) in enumerate(A):
        if k%10000 == 0:
            log.debug('Computing B^2_d. Element %s/%s.' % (k, len(A)))
        add_tup = True
        for t in range(d):
            if np.gcd(t,d) == 1:
                vec = [t*a/d for a in tup]
                vecNew = [entry-np.floor(entry) for entry in vec]
                if int(np.rint(sum(vecNew)))!=2:
                    add_tup = False
                    break
        if add_tup:
            B.add(tup)
    log.debug('len(B)=%s' % len(B))
    return B

def get_picard(matrix):
    '''
    Get the Picard number of a non-singular Delsarte K3 surface. Input is not
    checked for validity.
    :param matrix: A 4x4 matrix encoding a non-singular Delsarte K3 surface
    :return: The Picard number of the input matrix which is an integer between 1 and 20
    '''
    A = get_A_squared_d(matrix)
    B = get_B_squared_d(matrix)
    L = get_L_A_module(matrix)

    I = A - B
    lam = len(I.intersection(L))
    rho = 22 - lam
    return rho

def main(skip_to=0):
    '''
    Compute the Picard numbers of all smooth Delsarte surfaces and write the results
    to a log file. Can resume old computation by specifying at which index to start
    the computation.
    :param skip_to: Start computation at which index
    '''
    gen = get_iterator_of_A()
    for (k, A) in enumerate(gen):
        if k < skip_to:
            continue
        A = np.array(A)
        if check_valid_A(A):
            log.info('ID: %s, rho=%s, A=%s' % (k, get_picard(A), list(A)))


if __name__ == '__main__':
    main()
