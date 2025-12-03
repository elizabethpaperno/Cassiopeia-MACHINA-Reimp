import unittest
import pandas as pd
from cassiopeia_greedy import (
    compute_mutation_frequencies,
    best_mutation,
    split_on_mutation,
    greedy_reconstruct,
    SimpleNode,
)


class TestGreedySimplified(unittest.TestCase):

    def setUp(self):
        # Simple test matrix
        # Columns = sites, rows = samples (cell barcodes)
        self.matrix = pd.DataFrame(
            {
                0: [0, 1, 1, 0],
                1: [0, 0, 2, 2],
                2: [1, 1, 1, 1],
            },
            index=["A", "B", "C", "D"],
        )
        self.samples = ["A", "B", "C", "D"]
        self.missing = -1

    def test_compute_mutation_frequencies(self):
        freqs = compute_mutation_frequencies(self.matrix, self.samples, self.missing)
        # site 0 counts: 0→2, 1→2
        self.assertEqual(freqs[0][0], 2)
        self.assertEqual(freqs[0][1], 2)
        # site 1 counts: 0→2, 2→2
        self.assertEqual(freqs[1][2], 2)
        self.assertEqual(freqs[1][0], 2)
        # site 2 counts: only 1→4
        self.assertEqual(freqs[2][1], 4)

    def test_best_mutation(self):
        # site 2 should be ignored because mutation 1 occurs in ALL samples
        # site 0: mutation 1 occurs in 2 cells
        # site 1: mutation 2 occurs in 2 cells
        site, state = best_mutation(self.matrix, self.samples, self.missing)

        # Both sites have tied frequency, but site 0 < site 1 by index,
        # so site 0 wins because code loops sites in order.
        self.assertEqual(site, 0)
        self.assertEqual(state, 1)

    def test_best_mutation_with_weights(self):
        # Weight site 1 heavily
        weights = {
            0: {1: 1},
            1: {2: 10},
            2: {1: 1},
        }
        site, state = best_mutation(self.matrix, self.samples, self.missing, weights)

        # site 1 should now win because 2 cells × weight 10 = 20
        self.assertEqual(site, 1)
        self.assertEqual(state, 2)

    def test_split_on_mutation(self):
        # Split on site 0, state 1
        left, right = split_on_mutation(
            self.matrix, self.samples, site=0, state=1, missing_state_indicator=self.missing
        )
        self.assertCountEqual(left, ["B", "C"])
        self.assertCountEqual(right, ["A", "D"])

    def test_split_with_missing_values(self):
        # Inject missing data
        mat2 = self.matrix.copy()
        mat2.loc["A", 0] = -1
        left, right = split_on_mutation(
            mat2, self.samples, 0, 1, missing_state_indicator=-1
        )
        # Missing sample A joins the bigger group
        # site 0 values: A=?, B=1, C=1, D=0  ⇒ left has 2, right has 1
        self.assertIn("A", left)

    def test_greedy_reconstruct_structure(self):
        tree = greedy_reconstruct(self.matrix, self.samples)

        # Root should split into two children
        self.assertFalse(tree.is_leaf())
        self.assertEqual(len(tree.children), 2)

        # Leaves should be SimpleNodes with actual sample lists
        leaves = []

        def collect_leaves(node):
            if node.is_leaf():
                leaves.append(node.samples)
            else:
                for child in node.children:
                    collect_leaves(child)

        collect_leaves(tree)

        # Whatever the exact splits, all samples must appear exactly once total
        flat = [s for group in leaves for s in group]
        self.assertCountEqual(flat, self.samples)

    def test_leaf_when_no_split(self):
        # All samples have identical states → no valid mutation → leaf
        mat = pd.DataFrame(
            {
                0: [0, 0, 0],
                1: [1, 1, 1],
            },
            index=["X", "Y", "Z"],
        )
        tree = greedy_reconstruct(mat, ["X", "Y", "Z"])

        self.assertTrue(tree.is_leaf())
        self.assertCountEqual(tree.samples, ["X", "Y", "Z"])


if __name__ == "__main__":
    unittest.main()
