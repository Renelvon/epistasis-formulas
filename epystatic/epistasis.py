"""
COMPUTE EPISTASIS
=================
"""

import itertools

import numpy as np

from epystatic import circuits
from epystatic import fourier
from epystatic import slicing
from epystatic import utils


def compute_epistasis(w, w_err, species=5, interval=False):
    """Compute all epistatic interactions, via circuits and interaction
    coordinates among 'species' species. Return calculations in the form of
    3 matrices: results, standard errors, calculation tag.

    If 'interval' is False, then results and standard errors are returned;
    otherwise, low and high limits of epistasis are returned (intervals).
    """
    # First compute all circuit interactions
    c_epi = gen_all_circuit_interactions(w, w_err, species, interval=interval)
    c_results, c_errors, c_tags = c_epi
    results, errors, tags = [c_results], [c_errors], [c_tags]

    # Second, compute all significant orders of interaction coordinates
    u_epi = gen_all_coordinate_interactions(
        w, w_err, species, interval=interval
    )
    for u_results, u_errors, u_tags in u_epi:
        results.append(u_results)
        errors.append(u_errors)
        tags.append(u_tags)

    # Finally, join results and return to user
    all_results = np.hstack(results)
    all_errors = np.hstack(errors)
    all_tags = tuple(itertools.chain.from_iterable(tags))
    return all_results, all_errors, all_tags


def gen_all_circuit_interactions(w, w_err, species, interval=False):
    """Compute all circuit interactions and standard errors. Use interval
    arithmetic if interval==True.

    If 'interval' is True, then w is interpreted as the fitness vector
    corresponding to the low end of the intervals and w_err as the fitness
    vector corresponding to the high end of intervals. Instead of results and
    errors, intervals of epistasis are computed.
    """
    assert 2 ** species == len(w), 'Fitness vector has bad length'
    assert 2 ** species == len(w_err), 'Fitness error vector has bad length'

    # Generate all circuit interactions of order 3
    # Shape: (num_of_circtuits, 2**3)
    circuits_a2t = circuits.gen_circuits_3()

    # Generate all possible setups (projections) of order 3
    # Shape: (num_of_projections, 2**3)
    wp, w_tags = slicing.generate_all_projections(w, 2, species, 3)
    wp_err = slicing.generate_all_projections(
        w_err, 2, species, 3, tagged=False
    )

    if interval:
        w_low_p = w
        w_high_p = w_err

        # Split circuit matrix into nonnegative and nonpositive part,
        # and use elementary rules of interval arithmetic to calcualte results.
        pos_a2t, neg_a2t = utils.split_pos_neg(circuits_a2t)

        results_low = pos_a2t.dot(w_low_p) + neg_a2t.dot(w_high_p)
        results_high = pos_a2t.dot(w_high_p) + neg_a2t.dot(w_low_p)

        results = results_low
        errors = results_high
    else:
        # Calculate all interactions via inner product of circuits and
        # fitnesses (wp).
        # Result shape: (num_of_circuits, num_of_projections)
        results = circuits_a2t.dot(wp.T)

        # Calculate std error of all interactions via inner product of the
        # squares of coefficients (circuit matrix) and squares of fitness
        # errors (wp_err).
        # Result shape: (num_of_circuits, num_of_projections)
        errors = np.sqrt((circuits_a2t ** 2).dot((wp_err ** 2).T))

    # Generate tags for all circuit interactions.
    all_c_tags = (
        gen_circuit_tag(i, w_tag)
        for i in range(len(circuits_a2t))
        for w_tag in w_tags
    )

    # Remove duplicate interactions from the results (according to the
    # interaction formula computed; see 'circuits').
    return filter_duplicates(results, errors, all_c_tags)


def filter_duplicates(results, errors, tags):
    """Drop results and corresponding errors if they correspond to known
    duplicate circuit calculations.

    For more information, see 'circuits'.
    """
    unique_tags, drop_indices = [], []
    for i, tag in enumerate(tags):
        if tag in circuits.DUPLICATES:
            drop_indices.append(i)
        else:
            unique_tags.append(tag)

    unique_results = np.delete(results.flatten(), drop_indices)
    unique_errors = np.delete(errors.flatten(), drop_indices)

    return unique_results, unique_errors, unique_tags


def gen_circuit_tag(i, w_tag):
    """Generate a tag corresponding to a circuit interaction.

    It is the concatenation of the circuit name along with the projection tag.
    """
    return '{0:c}_{1:s}'.format(ord('a') + i, w_tag)


def gen_all_coordinate_interactions(w, w_err, species, interval=False):
    """Generate all non-singleton coordinate interactions and standard errors
    of all significant orders 'order' among 'species' species.

    If 'interval' is True, then w is interpreted as the fitness vector
    corresponding to the low end of the intervals and w_err as the fitness
    vector corresponding to the high end of intervals. Instead of results and
    errors, intervals of epistasis are computed.
    """
    assert 2 ** species == len(w), 'Bad fitness vector length'
    assert 2 ** species == len(w_err), 'Bad fitness std vector length'

    for order in range(3, species + 1):  # order==2 is covered by circuits
        yield gen_coordinate_interactions_of_order(
            w, w_err, species, order, interval=interval
        )


def gen_coordinate_interactions_of_order(
    w, w_err, species, order, interval=False
):
    """Generate all non-singleton coordinate interactions and standard errors
    of given 'order' among 'species' species.

    If 'interval' is True, then w is interpreted as the fitness vector
    corresponding to the low end of the intervals and w_err as the fitness
    vector corresponding to the high end of intervals. Instead of results and
    errors, intervals of epistasis are computed.
    """
    # Shape: (2**order - order - 1, 2**order)
    f_mat = fourier.generate_fourier_matrix(order)

    # Generate all possible setups (projections) of given order
    # Shape: (num_of_projections, 2**order)
    wp, w_tags = slicing.generate_all_projections(w, 2, species, order)
    wp_err = slicing.generate_all_projections(
        w_err, 2, species, order, tagged=False
    )

    if interval:
        w_low_p = w
        w_high_p = w_err

        # Split circuit matrix into nonnegative and nonpositive part,
        # and use elementary rules of interval arithmetic to calcualte results.
        pos_a2t, neg_a2t = utils.split_pos_neg(f_mat)

        results_low = pos_a2t.dot(w_low_p) + neg_a2t.dot(w_high_p)
        results_high = pos_a2t.dot(w_high_p) + neg_a2t.dot(w_low_p)

        results = results_low
        errors = results_high
    else:
        # Calculate all interactions via inner product
        # Shape: (2**order - order - 1, num_of_projections)
        results = f_mat.dot(wp.T)

        # Calculate std error of all interactions
        # Shape: (2**order - order - 1, num_of_projections)
        errors = np.sqrt((f_mat ** 2).dot((wp_err ** 2).T))

    u_tags = tuple(
        gen_coordinate_tag(i, order, w_tag)
        for i in range(2 ** order)
        for w_tag in w_tags
        if i & (i - 1)  # evaluates to 0 iff i == 0 or i == 2**j for j >= 0
    )

    return results.flatten(), errors.flatten(), u_tags


# TODO: Write more clearly
# TODO: Optimize
def gen_coordinate_tag(formula, order, w_tag):
    """Generate a tag corresponding to an interaction coordinate

    The interpretation of the tag is as follows: Numbers appear in the
    position of fixed (background) species, denoting the fixed value; letters
    appear in the non-fixed species, in capital or small form, depending on the
    index of the row of the Fourier matrix that is applied to the projection.

    For example, the interaction u_101 applied to the w_A1BC0 projection will
    result in the u_A1bC0 result tag, as A and C correspond to '1' in u_101,
    whereas B corresponds to 0.
    """
    aux = 1 << order
    desc = list(w_tag)
    for i, c in enumerate(desc):
        if c not in '01':
            aux >>= 1
            if not formula & aux:
                desc[i] = c.lower()
    return 'u_' + ''.join(desc)


if __name__ == '__main__':
    # w_data = np.array([32] + list(range(1, 32)))
    w_data = np.array(
        [
            10.455,
            11.167,
            10.125,
            10.542,
            10.458,
            9.875,
            10.167,
            10.542,
            10.375,
            9.833,
            9.875,
            10.292,
            9.625,
            10.125,
            9.958,
            9.792,
            9.833,
            10.167,
            9.875,
            10.250,
            10.125,
            9.833,
            9.917,
            9.583,
            9.958,
            10.042,
            10.000,
            9.708,
            9.875,
            10.042,
            9.792,
            10.125,
        ]
    )
    w_err = np.array(
        [
            0.143,
            0.177,
            0.110,
            0.134,
            0.104,
            0.125,
            0.115,
            0.134,
            0.132,
            0.098,
            0.092,
            0.127,
            0.118,
            0.125,
            0.165,
            0.134,
            0.130,
            0.167,
            0.110,
            0.138,
            0.125,
            0.130,
            0.133,
            0.119,
            0.127,
            0.112,
            0.147,
            0.153,
            0.125,
            0.127,
            0.104,
            0.110,
        ]
    )

    w_pystasis = utils.convert_vector_to_pystasis_order(w_data)
    w_pystasis_err = utils.convert_vector_to_pystasis_order(w_err)
    w = w_pystasis
    w_err = w_pystasis_err
    results, errors, tags = compute_epistasis(w, w_err, 5)
    for result, error, tag in zip(results, errors, tags):
        print(tag, result, error)
