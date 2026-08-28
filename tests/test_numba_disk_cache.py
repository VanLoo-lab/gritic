import unittest

from numba.core.registry import CPUDispatcher

from gritic import hitandrun, sampletools


class NumbaDiskCacheTest(unittest.TestCase):
    def test_every_numba_kernel_uses_the_disk_cache(self):
        kernels = []
        for module in (hitandrun, sampletools):
            kernels.extend(
                value
                for value in vars(module).values()
                if isinstance(value, CPUDispatcher)
            )

        self.assertGreater(len(kernels), 0)
        for kernel in kernels:
            with self.subTest(kernel=kernel.py_func.__name__):
                self.assertTrue(kernel._cache._enabled)


if __name__ == '__main__':
    unittest.main()
