"""
CIRCUIT GENERATION
==================

This module generates the circuit matrix C of all the circuits from a to t,
as defined in [1]. No matter to what data they are applied, the circuits
calculated here always calcuate interactions between 3 active species, and thus
are linear combinations of fitness vectors of length 2^3==8. Thus, applying
them on data of more than 3 species requires us to select the active (non-fixed)
species and the background (fixed) species, project the data and calculate the
resulting epistases.

Although all circuits are independent when applied to data of 3 species, this
is not the case when calculating circuits over projected higher-order data. A
suitable function is provided that generates duplicate circuits for the hard-
coded case of calculating all order-3 circuits when n==5 species are available.

To avoid extraneous computations, the duplicate list is provided.

[1] Epistasis and Shapes of Fitness Landscapes
"""

import itertools

import numpy as np

from epystatic import fourier
from epystatic import slicing
from epystatic import utils

DUPLICATES = set(
    (
        'a_0ABC0',
        'a_1ABC0',
        'a_A0BC0',
        'a_A1BC0',
        'a_AB0C0',
        'a_AB1C0',
        'a_ABC00',
        'a_ABC01',
        'a_ABC10',
        'b_0AB0C',
        'b_0ABC0',
        'b_0ABC1',
        'b_1AB0C',
        'b_1ABC0',
        'b_1ABC1',
        'b_A0B0C',
        'b_A0BC0',
        'b_A0BC1',
        'b_A1B0C',
        'b_A1BC0',
        'b_A1BC1',
        'b_AB00C',
        'b_AB01C',
        'b_AB0C0',
        'b_AB0C1',
        'b_AB10C',
        'b_AB1C0',
        'b_AB1C1',
        'b_ABC00',
        'b_ABC01',
        'b_ABC10',
        'b_ABC11',
        'c_0AB0C',
        'c_0ABC0',
        'c_0ABC1',
        'c_1AB0C',
        'c_1ABC0',
        'c_1ABC1',
        'c_A0B0C',
        'c_A0BC0',
        'c_A0BC1',
        'c_A1B0C',
        'c_A1BC0',
        'c_A1BC1',
        'c_AB00C',
        'c_AB01C',
        'c_AB0C0',
        'c_AB0C1',
        'c_AB10C',
        'c_AB1C0',
        'c_AB1C1',
        'c_ABC00',
        'c_ABC01',
        'c_ABC10',
        'c_ABC11',
        'd_0A0BC',
        'd_0AB0C',
        'd_0AB1C',
        'd_0ABC0',
        'd_0ABC1',
        'd_1A0BC',
        'd_1AB0C',
        'd_1AB1C',
        'd_1ABC0',
        'd_1ABC1',
        'd_A00BC',
        'd_A01BC',
        'd_A0B0C',
        'd_A0B1C',
        'd_A0BC0',
        'd_A0BC1',
        'd_A10BC',
        'd_A1B0C',
        'd_A1B1C',
        'd_A1BC0',
        'd_A1BC1',
        'd_AB00C',
        'd_AB01C',
        'd_AB0C0',
        'd_AB0C1',
        'd_AB10C',
        'd_AB11C',
        'd_AB1C0',
        'd_AB1C1',
        'd_ABC00',
        'd_ABC01',
        'd_ABC10',
        'd_ABC11',
        'e_0A0BC',
        'e_0AB0C',
        'e_0AB1C',
        'e_0ABC0',
        'e_0ABC1',
        'e_1A0BC',
        'e_1AB0C',
        'e_1AB1C',
        'e_1ABC0',
        'e_1ABC1',
        'e_A00BC',
        'e_A01BC',
        'e_A0B0C',
        'e_A0B1C',
        'e_A0BC0',
        'e_A0BC1',
        'e_A10BC',
        'e_A1B0C',
        'e_A1B1C',
        'e_A1BC0',
        'e_A1BC1',
        'e_AB00C',
        'e_AB01C',
        'e_AB0C0',
        'e_AB0C1',
        'e_AB10C',
        'e_AB11C',
        'e_AB1C0',
        'e_AB1C1',
        'e_ABC00',
        'e_ABC01',
        'e_ABC10',
        'e_ABC11',
        'f_00ABC',
        'f_01ABC',
        'f_0A0BC',
        'f_0A1BC',
        'f_0AB0C',
        'f_0AB1C',
        'f_0ABC0',
        'f_0ABC1',
        'f_10ABC',
        'f_1A0BC',
        'f_1A1BC',
        'f_1AB0C',
        'f_1AB1C',
        'f_1ABC0',
        'f_1ABC1',
        'f_A00BC',
        'f_A01BC',
        'f_A0B0C',
        'f_A0B1C',
        'f_A0BC0',
        'f_A0BC1',
        'f_A10BC',
        'f_A11BC',
        'f_A1B0C',
        'f_A1B1C',
        'f_A1BC0',
        'f_A1BC1',
        'f_AB00C',
        'f_AB01C',
        'f_AB0C0',
        'f_AB0C1',
        'f_AB10C',
        'f_AB11C',
        'f_AB1C0',
        'f_AB1C1',
        'f_ABC00',
        'f_ABC01',
        'f_ABC10',
        'f_ABC11',
    )
)


def gen_circuits_3():
    """Generate the circuit matrix C.

    It is a 20 x 8 matrix which is derived from the full Fourier interaction
    matrix of order 3, as specified in page 1325 of [1]. Note that computing
    the circuit matrix via the Fourier matrix is a matter of convenience, but
    causes a lot of numerical cancellations to happen. Thus, to avoid subtle
    numerical errors and loss of precision in the final results, C should be
    computed *prior* to applying it to data.
    """
    f3 = fourier.generate_fourier_matrix(3)
    return np.vstack(
        (gen_circuits_3_a2f(f3), gen_circuits_3_g2l(f3), gen_circuits_3_m2t(f3))
    )


def gen_circuits_3_a2f(f_mat):
    """Generate circuits a to f, using the formulas of [1]."""
    u_full = f_mat[3]  # select 'u_111'
    circuits = np.empty((6, 8), dtype=f_mat.dtype)
    idx = 0
    for u_row in reversed(f_mat[:3]):  # for each of [u_011, u_101, u_110]...
        circuits[idx] = u_row + u_full
        circuits[idx + 1] = u_row - u_full
        idx += 2
    return circuits // 2  # compensate


def gen_circuits_3_g2l(f_mat):
    """Generate circuits g to l, using the formulas of [1]."""
    circuits = np.empty((6, 8), dtype=f_mat.dtype)
    combs = list(itertools.combinations(range(3), 2))
    combs.reverse()
    idx = 0

    # for each ordered pair of [u_011, u_101, u_110]...
    for u_low, u_high in combs:
        circuits[idx] = f_mat[u_high] + f_mat[u_low]
        circuits[idx + 1] = f_mat[u_high] - f_mat[u_low]
        idx += 2
    return circuits // 2  # compensate


def gen_circuits_3_m2t(f_mat):
    """Generate circuits m to t, using the formulas of [1]."""
    sign_m = np.array(
        [
            [-1, -1, -1, -1],
            [-1, -1, -1, 1],
            [1, 1, -1, -1],
            [1, 1, -1, 1],
            [1, -1, 1, -1],
            [1, -1, 1, 1],
            [-1, 1, 1, 1],
            [-1, 1, 1, -1],
        ]
    )
    # Combine elements of the Fourier matrix to form last 8 circuits
    # (caution: cancellations!)
    return sign_m.dot(f_mat) // 2  # compensate


def generate_duplicates():
    """Generate circuits m to t, using the formulas of [1]."""
    c = gen_circuits_3()

    # Create a special fitness vector of size 32, where w_j == j.
    w = np.arange(0, 2 ** 5)

    # Generate all tagged projections of order 3 from 5 species.
    wp, wtags = slicing.generate_all_projections(w, 2, 5, 3)

    # Calculate the (unreduced) result of applying the circuit matrix
    # to each projection. Because of the specific W chosen, this
    # results in vectors that contain signed elements of W and 0s.
    epistasis_setups = (
        (c[i] * wp[j], utils.gen_tag(i, wtags[j]))
        for i, j in itertools.product(range(len(c)), range(len(wp)))
    )

    # Drop 0s from all resulting vectors, convert them to tuples and produce a
    # dictionary as follows: The first time an epistatic calculation produces a
    # tuple, add the tuple to the keys, without storing the epistatic tag. For
    # each subsequent tag which recalcualtes the same tuple, record the tag.
    # Finally, output all values as computing duplicates.
    #
    # NOTE: We assume that the tuples produced are pre-sorted in ascending
    # magnitude of elements. In our implementation this is guaranteed by using
    # the specific W above and the NumPy projecting facilities via the
    # 'slicing' module.
    duplicates = {}
    for e, ctx in epistasis_setups:
        enz = tuple(e[e.nonzero()])
        if enz in duplicates:
            duplicates[enz].append(ctx)
        else:
            duplicates[enz] = []

    return set(itertools.chain.from_iterable(duplicates.values()))


if __name__ == '__main__':
    print('Circuits')
    c_idx = ord('a')
    for c in gen_circuits_3():
        print(chr(c_idx), c)
        c_idx += 1

    print('\nDuplicates when applying circuits to 5-to-3 binary projections')
    print(sorted(generate_duplicates()))
