import unittest

import numpy as np

from epystatic import fourier


class TestFourierMatrixGeneration(unittest.TestCase):

    @staticmethod
    def _inner_product(i, j):
        return sum(1 for s in np.binary_repr(i & j) if s == '1')

    @classmethod
    def _generate_full_fourier_matrix_iter(cls, n):
        size = 2**n
        f = np.empty((size, size), dtype=np.int16)
        for i in range(size):
            for j in range(size):
                if cls._inner_product(i, j) % 2 == 0:
                    f[i, j] = 1
                else:
                    f[i, j] = -1
        return f

    def test_generate_full_fourier_matrix_0(self):
        self.assertTrue(np.array_equal(
            fourier.generate_full_fourier_matrix(0),
            self._generate_full_fourier_matrix_iter(0),
        ))

    def test_generate_full_fourier_matrix_1(self):
        self.assertTrue(np.array_equal(
            fourier.generate_full_fourier_matrix(1),
            self._generate_full_fourier_matrix_iter(1),
        ))

    def test_generate_full_fourier_matrix_2(self):
        self.assertTrue(np.array_equal(
            fourier.generate_full_fourier_matrix(2),
            self._generate_full_fourier_matrix_iter(2),
        ))

    def test_generate_full_fourier_matrix_5(self):
        self.assertTrue(np.array_equal(
            fourier.generate_full_fourier_matrix(5),
            self._generate_full_fourier_matrix_iter(5),
        ))

    def test_generate_singleton_indices_0(self):
        self.assertCountEqual(
            fourier.generate_singleton_indices(0),
            (0,)
        )

    def test_generate_singleton_indices_1(self):
        self.assertCountEqual(
            fourier.generate_singleton_indices(1),
            (0,1)
        )

    def test_generate_singleton_indices_5(self):
        self.assertCountEqual(
            fourier.generate_singleton_indices(5),
            (0, 1, 2, 4, 8, 16)
        )

    def test_generate_fourier_matrix_2(self):
        self.assertTrue(np.array_equal(
            fourier.generate_fourier_matrix(2),
            np.array([[1, -1, -1, 1]])
        ))

    def test_generate_fourier_matrix_3(self):
        self.assertTrue(np.array_equal(
            fourier.generate_fourier_matrix(3),
            np.array([
                [ 1, -1, -1,  1,  1, -1, -1,  1],
                [ 1, -1,  1, -1, -1,  1, -1,  1],
                [ 1,  1, -1, -1, -1, -1,  1,  1],
                [ 1, -1, -1,  1, -1,  1,  1, -1]
            ])
        ))


if __name__ == '__main__':
    unittest.main()
