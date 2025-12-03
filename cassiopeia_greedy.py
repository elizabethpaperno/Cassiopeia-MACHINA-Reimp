import numpy as np
import pandas as pd
from collections import defaultdict


class SimpleNode:
    """Node for tree"""
    def __init__(self, name=None, samples=None):
        self.name = name
        self.children = []
        self.samples = samples if samples is not None else []

    def add_child(self, child):
        self.children.append(child)

    def is_leaf(self):
        return len(self.children) == 0


def compute_mutation_frequencies(character_matrix, samples, missing_state_indicator):
    """
    Computes how frequently each state occurs at each character (column).
    Returns: {col: {state: count}}
    """
    df = character_matrix.loc[samples]
    num_sites = df.shape[1]

    freqs = {col: defaultdict(int) for col in range(num_sites)}
    for _, row in df.iterrows():
        for col, state in enumerate(row):
            freqs[col][state] += 1

    return freqs


def best_mutation(character_matrix, samples, missing_state_indicator, weights=None):
    """
    Select (site, state) with highest frequency (optionally weighted).
    Mimics Cassiopeia's logic: choose the mutation that occurs most often,
    excluding 0 and missing-state.
    """
    frequencies = compute_mutation_frequencies(
        character_matrix, samples, missing_state_indicator
    )

    best_score = 0
    chosen_site = None
    chosen_state = None

    for site in frequencies:
        for state, count in frequencies[site].items():
            if state == missing_state_indicator or state == 0:
                continue  # skip missing and reference state

            # Avoid splitting on mutations shared by ALL samples
            if count == len(samples):
                continue

            score = count * (weights[site][state] if weights else 1)

            if score > best_score:
                best_score = score
                chosen_site = site
                chosen_state = state

    return chosen_site, chosen_state


def split_on_mutation(character_matrix, samples, site, state, missing_state_indicator):
    """
    Splits samples into left (have mutation) and right (do not).
    Missing data go into their own buffer.
    No fancy missing data resolution — simplified logic.
    """
    left = []
    right = []
    missing = []

    df = character_matrix.loc[samples]

    for sample in samples:
        value = df.loc[sample, site]
        if value == missing_state_indicator:
            missing.append(sample)
        elif value == state:
            left.append(sample)
        else:
            right.append(sample)

    # SIMPLE strategy: assign missing to the larger side
    if missing:
        if len(left) >= len(right):
            left.extend(missing)
        else:
            right.extend(missing)

    return left, right


def greedy_reconstruct(
    character_matrix,
    samples=None,
    name="root",
    missing_state_indicator=-1,
    weights=None,
):
    """
    Main tree reconstruction recursion, using Cassiopeia-style inputs.

    character_matrix: DataFrame of mutation states.
    samples: list of sample names to process.
    """
    if samples is None:
        samples = list(character_matrix.index)

    # Base case
    if len(samples) <= 1:
        return SimpleNode(name=name, samples=samples)

    # Choose split mutation
    site, state = best_mutation(
        character_matrix, samples, missing_state_indicator, weights
    )

    # If no valid split: leaf
    if site is None:
        return SimpleNode(name=name, samples=samples)

    # Partition samples
    left, right = split_on_mutation(
        character_matrix, samples, site, state, missing_state_indicator
    )

    if len(left) == 0 or len(right) == 0:
        return SimpleNode(name=name, samples=samples)

    # Make the parent node
    parent = SimpleNode(name=f"{name}_site{site}_state{state}", samples=[])

    left_child = greedy_reconstruct(
        character_matrix,
        samples=left,
        name=f"{name}_L",
        missing_state_indicator=missing_state_indicator,
        weights=weights,
    )
    right_child = greedy_reconstruct(
        character_matrix,
        samples=right,
        name=f"{name}_R",
        missing_state_indicator=missing_state_indicator,
        weights=weights,
    )

    parent.add_child(left_child)
    parent.add_child(right_child)
    return parent


def print_tree(node, indent=0):
    """Print tree"""
    prefix = " " * indent
    if node.is_leaf():
        print(f"{prefix}- Leaf: samples={node.samples}")
    else:
        print(f"{prefix}- Node: {node.name}")
        for child in node.children:
            print_tree(child, indent + 4)
