"""
DATA RESHAPING AND PROJECTING (SLICING)
=======================================

This module contains functions that reshape the input data (fitness vector W)
into a high-dimensional tensor. This allows to utilize the facilities of NumPy
to systematically select the proper elements to use when calculating epistatic
effects against various backgrounds.

The basis of this module's functionality is the following: assume that we are
given a fitness vector W of length 2^n (i.e. involving n species). The indices
of this vector denote setups, as explained in the 'fourier' module. To make the
interpretation of the w_j more explicit, we can reshape W into a tensor with n
axes, were each axes has dimension 2, that is, addressed by a single bit. This
shape allows us to group the fitness values in a structure that ---albeit more
complicated at first glance--- is actually much easier to use via NumPy.

As an example of what we intend to do, consider the following case: We are
given a 1-D vector W of length 32, corresponding to fitnesses of 5-species
setups. Suppose we want to access the fitness of the setup where species 2, 4,
5 are present, and 1 and 3 are absent. This corresponds to the binary string
11010, that is, number 26. We can access this value via

    W[26]

Let us transform W into a 5-axes tensor, V. The required value can now be
accessed via

    V[1][1][0][1][0]

that is, by directly using the setup's bitstring to access the tensor, without
any calculations.

The true power of the tensor representation comes into play when we are trying
to calculate low-order effects (e.g. order-3 circuits) against a specified
background. We can then use the excellent NumPy slicing facilities to
automatically select the right elements out of V. For example: Among 3 species,
the value of circuit A is the following:

    a = w_000 - w_010 - w_100 + w_110

This implies a vector W of length 8. For 4 indices of W the cofactor is 0, for
2 more elements the cofactor is +1 and for the rest it is -1. Suppose now W is
of length 32 (all 5 species) and we want to calculate a similar circuit among
species 2, 4 and 5 when species 1 and 3 form the background, where species 1
is absent and species 3 is present. This would require us to calculate the
formula above but with values of the form

    w_XY1Z0

to correctly incorporate the background. Selecting the proper 8 values out of W
is non-trivial and slow in the 1D representation, as we would have to
interconvert between binary strings and decimals all the time. The same job in
the tensor representation becomes much easier, because NumPy allows us to
'slice' (project) the 5D tensor to a proper 3D tensor, which contains all the
elements we need, and then extract these elements in proper order. Essentially,
we form the expression

    V[*][*][1][*][0]

where '*' means 'all elements along this axis'. Once calculated, this results
in a 3D tensor, V_3. We then proceed to reshape V_3 into a 1D vector W_3 which
contains exactly the values W_XY1Z0, where XYZ are replaced by all 8 elements
of {0, 1}^3, in the correct order. More details below.
"""

import itertools

import numpy as np

from epystatic import utils


class TensorProjector:
    """
    An object that allows manipulation and slicing of a 1D vector as
    a multidimensional tensor.
    """

    def __init__(self, base, rank):
        """
        Creates an object that knows how to manipulate a 1D vector of
        length base^rank, corresponding to the respective combinatorial setups.

        We usually expect base==2, as there are 2 options per species in our
        experiment. Similarly, we expect rank==5, as this is the number of
        species. None of these assumptions is strict in this module.
        """
        self.base = base
        self.rank = rank

        """This is the '*' element, which can be used to select all dimensions
        along an axis (all options for a given species)."""
        self.all = slice(0, base, 1)

    def tensorize(self, v):
        # TODO: Use values.reshape for panda series
        return v.reshape(*[self.base] * self.rank)

    def project_vector(self, v, const_indices):
        """Select the appropriate elements from a 1D vector 'v', as
        corresponding to the background specified by 'const_indices'.

        'const_indices' should be a sequence of 2-tuples of the form (i, q)
        where the value q is specified for index i. This is interpreted as
        selecting fitnesses where species i always has value q. Multiple
        fixed axes are allowed.
        """
        v_tensor = self.tensorize(v)
        return self.project_tensor(v_tensor, const_indices)

    def project_tensor(self, v, const_indices):
        projection_desc = [self.all] * self.rank
        for index, value in const_indices:
            projection_desc[index] = value
        return v[projection_desc]


def generate_all_standard_projections(
    w, base, full_rank, proj_rank, tagged=True
):
    """Systematically generate all low-order 'projections' of a fitness vector
    W, against 'standard' (i.e. null, all-0) backgrounds. Generate tags when
    'tagged' is True.

    This function accepts a 1D fitness vector W corresponding to all possible
    experimental setups of 'full_rank' species, where each species is assumed
    to have 'base' options. It also accepts a 'proj_rank' number, which is the
    number of species that are considered non-fixed. According to the above,
    it generates all possible 'contexts' (backgrounds) of fixed species,
    sets all fixed species to 0 (i.e. absent) and extracts from W the subset of
    relevant fitnesses, in the order they would appear if the experiment was
    carried out in such a reduced environment.

    It optionally tags the projections with names, so that they can be easily
    understood. For example, the tag 'A0B0C' corresponds to a (standard)
    projection where species 1, 3, 5 are non-fixed and 2 and 4 are fixed to 0
    (absent).
    """
    assert len(w) == base ** full_rank
    assert 0 < proj_rank <= full_rank

    res_rank = full_rank - proj_rank
    fixed_setup = (0,) * res_rank  # standard projections: absent bystanders
    idx_options = tuple(itertools.combinations(range(full_rank), r=res_rank))
    contexts = tuple(tuple(zip(io, fixed_setup)) for io in idx_options)

    tp = TensorProjector(base, full_rank)
    wt = tp.tensorize(w)

    # For each projection context, select and store the apppropriate elements
    # from fitness tensor
    projections = np.array(
        tuple(tp.project_tensor(wt, ctx).flatten() for ctx in contexts)
    )

    if not tagged:
        return projections

    # Create tags if requested. See above for interpretation of tag strings.
    # See 'utils' for explanation of context string.
    tags = tuple(utils.format_context(ctx, full_rank) for ctx in contexts)
    return projections, tags


def generate_all_projections(w, base, full_rank, proj_rank, tagged=True):
    """Systematically generate all low-order 'projections' of a fitness vector
    W, against all possible backgrounds. Generate tags when 'tagged' is True.

    This function accepts a 1D fitness vector W corresponding to all possible
    experimental setups of 'full_rank' species, where each species is assumed
    to have 'base' options. It also accepts a 'proj_rank' number, which is the
    number of species that are considered non-fixed. According to the above,
    it generates *all* possible 'contexts' (backgrounds) of fixed species,
    and extracts from W the subset of relevant fitnesses, in the order they
    would appear if the experiment was carried out in such a reduced
    environment.

    It optionally tags the projections with names, so that they can be easily
    understood. For example, the tag 'AB1C0' corresponds to a (non-standard)
    projection where species 2, 4, 5 are non-fixed, species 1 is fixed to 0 and
    species 3 is fixed to 1.
    """
    assert base > 0
    assert len(w) == base ** full_rank
    assert 0 < proj_rank <= full_rank

    res_rank = full_rank - proj_rank
    # (Non-)standard projections: all possible bystander backgrounds
    fixed_setups = tuple(itertools.product(range(base), repeat=res_rank))
    idx_options = tuple(itertools.combinations(range(full_rank), r=res_rank))
    contexts = tuple(
        tuple(zip(io, fixed_setup))
        for io, fixed_setup in itertools.product(idx_options, fixed_setups)
    )

    # Turn W into a tensor.
    tp = TensorProjector(base, full_rank)
    wt = tp.tensorize(w)

    # For each projection context, select and store the apppropriate elements
    # from fitness tensor
    projections = np.array(
        tuple(tp.project_tensor(wt, ctx).flatten() for ctx in contexts)
    )

    if not tagged:
        return projections

    # Create tags if requested. See above for interpretation of tag strings.
    tags = tuple(utils.format_context(ctx, full_rank) for ctx in contexts)
    return projections, tags


if __name__ == '__main__':
    print('Assume 5 different species, either present or absent (l = 2).')
    print('Construct the fitness vector w_flat, indexed by {b_0 b_1 ... b_4}')
    l, rank = 2, 5
    w_flat = np.arange(0, l ** rank)
    print(w_flat)

    print(
        '\nTransform w_flat into 5-D tensor with l = 2 elements per dimension.'
    )
    w_tensor = w_flat.reshape(*[l] * rank)
    print(w_tensor)

    print('\nSelect a 2x2 slice of w where (b_0, b_2, b_3) = (0, 1, 1)')
    w_proj = w_tensor[0, :, 1, 1, :]
    print(w_proj)

    print('\nSame procedure but using TensorProjector')
    const_dims = [(0, 0), (2, 1), (3, 1)]
    cp = TensorProjector(2, 5)
    w_tensor2 = cp.project_vector(w_flat, const_dims)
    print(w_tensor2)

    print('\nGenerate all standard projections of size 3)')
    projections, tags = generate_all_standard_projections(w_flat, l, rank, 3)
    for p, t in zip(projections, tags):
        print(t, p)

    print('\nGenerate all projections of size 3)')
    projections, tags = generate_all_projections(w_flat, l, rank, 3)
    for p, t in zip(projections, tags):
        print(t, p)
