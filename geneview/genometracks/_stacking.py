"""
Stacking algorithms for arranging overlapping genomic features.

Implements the same stacking modes as Gviz/UCSC Genome Browser:
    - full: Each feature gets its own row (no overlapping).
    - pack: Pack features tightly, but non-overlapping items share rows.
    - squish: Like pack but with reduced row height.
    - dense: All features collapsed into a single row (overlapping allowed).
    - hide: No features shown (returns all zeros).

These algorithms assign a stack row index (0-based) to each feature based
on its start/end coordinates, ensuring no two features in the same row overlap.
"""

import numpy as np
from typing import Optional


def compute_stacking(
    starts: np.ndarray,
    ends: np.ndarray,
    mode: str = "squish",
    min_distance: int = 1,
) -> np.ndarray:
    """Compute stack row assignments for a set of genomic features.

    Uses a greedy interval-scheduling algorithm: features are sorted by start
    position and assigned to the first available row where they don't overlap
    any existing feature.

    Parameters
    ----------
    starts : np.ndarray
        Array of feature start positions.
    ends : np.ndarray
        Array of feature end positions.
    mode : str
        Stacking mode: 'full', 'pack', 'squish', 'dense', or 'hide'.
    min_distance : int
        Minimum distance between features in the same row. Default is 1.

    Returns
    -------
    np.ndarray
        Array of stack row indices (0-based) for each feature.
    """
    n = len(starts)
    if n == 0:
        return np.array([], dtype=int)

    if mode == "hide":
        return np.zeros(n, dtype=int)

    if mode == "dense":
        return np.zeros(n, dtype=int)

    # For 'full', 'pack', 'squish' - use greedy interval scheduling
    return _greedy_stacking(starts, ends, min_distance)


def _greedy_stacking(
    starts: np.ndarray,
    ends: np.ndarray,
    min_distance: int = 1,
) -> np.ndarray:
    """Greedy interval scheduling for stack assignment.

    Sorts features by start position and assigns each to the first row
    where it doesn't overlap any previously assigned feature.

    Parameters
    ----------
    starts : np.ndarray
        Start positions.
    ends : np.ndarray
        End positions.
    min_distance : int
        Minimum gap between features in the same row.

    Returns
    -------
    np.ndarray
        Row assignments (0-based).
    """
    n = len(starts)
    stacks = np.zeros(n, dtype=int)

    # Sort indices by start position (then by end position for ties)
    order = np.lexsort((ends, starts))

    # Track the end position of the last feature in each row
    row_ends = []  # list of end positions for each row

    for idx in order:
        s = starts[idx]
        placed = False

        for row_idx, row_end in enumerate(row_ends):
            if s >= row_end + min_distance:
                # This feature fits in this row
                stacks[idx] = row_idx
                row_ends[row_idx] = ends[idx]
                placed = True
                break

        if not placed:
            # Need a new row
            stacks[idx] = len(row_ends)
            row_ends.append(ends[idx])

    return stacks


def get_num_stacks(stacks: np.ndarray) -> int:
    """Return the number of distinct stack rows used."""
    if len(stacks) == 0:
        return 0
    return int(stacks.max()) + 1


def get_stack_heights(
    stacks: np.ndarray,
    mode: str = "squish",
    stack_height_frac: float = 0.75,
    reverse_stacking: bool = False,
) -> dict:
    """Compute the vertical position and height for each stack row.

    Parameters
    ----------
    stacks : np.ndarray
        Stack row assignments from compute_stacking.
    mode : str
        Stacking mode.
    stack_height_frac : float
        Fraction of available height each row occupies.
    reverse_stacking : bool
        If True, reverse the y-order of stack rows so that the first row
        is at the bottom instead of the top. Default is False.

    Returns
    -------
    dict
        Keys are 'y_positions' (array of y-center for each row),
        'row_height' (height of each row in data coordinates),
        'total_rows' (number of rows used).
    """
    total_rows = get_num_stacks(stacks)
    if total_rows == 0:
        return {"y_positions": np.array([]), "row_height": 0, "total_rows": 0}

    if mode == "dense":
        # All in one row
        return {
            "y_positions": np.array([0.5]),
            "row_height": 1.0 * stack_height_frac,
            "total_rows": 1,
        }

    row_height = stack_height_frac / total_rows
    y_positions = np.linspace(
        1.0 - row_height / 2,
        row_height / 2,
        total_rows,
    )

    if reverse_stacking:
        y_positions = y_positions[::-1]

    return {
        "y_positions": y_positions,
        "row_height": row_height,
        "total_rows": total_rows,
    }
