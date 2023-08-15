import functools
from math import log
from interval import interval
from .isotopes import ISOTOPES

from tqdm.auto import tqdm

# Monoisotopic masses of common adduct ions
COMMON_PRECURSORS = [
    # Source: https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
    0.0,  # Nothing
    1.007825,  # H+
]

# Minimum MW of what can be considered a fragment
MIN_CHUNK = 5.0


def upper_bound_by_MW(mw):
    return max(0.075 * mw - 1, 0.0)


def lower_bound_by_MW(mw):
    return log(max(mw * 0.075, 1.0), 2)


def unify_trees(trees: list[dict]):
    """
    Recursively merge `trees` into one tree.
    """
    if not trees:
        return {}
    elif len(trees) == 1:
        return trees[0]
    else:
        child1, child2, *rest = trees
        child1_keys = set(child1 or {})
        child2_keys = set(child2 or {})
        common_keys = child1_keys.intersection(child2_keys)
        return {
            **{k: child1[k] for k in child1_keys - common_keys},
            **{k: child2[k] for k in child2_keys - common_keys},
            **{k: unify_trees([child1[k], child2[k]]) for k in common_keys},
        }


class MAEstimator:
    def __init__(self, same_level=True, tol=0.01, adduct_masses=COMMON_PRECURSORS):
        self.same_level = same_level
        self.tol = tol
        self.adduct_masses = adduct_masses

    @functools.cache
    def estimate_by_MW(self, mw):
        lower, upper = mw - self.tol, mw + self.tol
        for isotope, weight in ISOTOPES.items():
            if lower < weight < upper:
                # MW matches an isotope; MA = 0
                print(f"HIT: {mw} ~ {isotope} ({weight})")
                return interval([0.0, 0.0])
        return interval([lower_bound_by_MW(mw), upper_bound_by_MW(mw)])

    def estimate_MA(self, tree: dict[float, dict], mw: float, progress=False):
        children = unify_trees([tree.get(mw, None) or self.precursors(tree, mw)])
        child_estimates = {mw: self.estimate_by_MW(mw)}

        for child in tqdm(children) if progress else children:
            complement = mw - child
            if complement < MIN_CHUNK or child < MIN_CHUNK:
                continue

            common = [
                p
                for p in self.common_precursors(children, child, complement)
                if p > MIN_CHUNK
            ]

            if common:
                print(f"Common precursors of {child} & {complement}: {common}")

            # Simple child + complement with no common precursors
            ma_candidates = [
                self.estimate_MA(children, child)
                + self.estimate_MA(children, complement)
                + 1
            ]

            for precursor in common:
                chunks = [child - precursor, complement - precursor, precursor]
                if min(chunks) < MIN_CHUNK:
                    continue
                chunk_mas = sum(
                    self.estimate_MA(
                        children,
                        chunk,
                    )
                    for chunk in chunks
                )
                ma_candidates.append(chunk_mas + 3)

            child_estimates[child] = min(ma_candidates, key=lambda x: x[0].sup)

        estimate = min(child_estimates.values(), key=lambda x: x[0].sup)
        return estimate

    def common_precursors(self, data, parent1, parent2):
        precursors1 = self.precursors(data, parent1)
        precursors2 = self.precursors(data, parent2)
        return set(precursors1).intersection(precursors2)

    def precursors(self, data, parent):
        if parent < MIN_CHUNK:
            return {}

        possible_ions = [parent + adduct for adduct in self.adduct_masses]
        parent_candidates = [
            d
            for d in data
            if any(d - self.tol < p < d + self.tol for p in possible_ions)
        ]
        children = unify_trees(
            [
                {
                    **{
                        p - child: self.same_level_precursors(data, p - child)
                        for child in data[p] or {}
                    },
                    **(data[p] or {}),
                }
                for p in parent_candidates
            ]
        )
        if not children and self.same_level:
            children = unify_trees(
                [
                    self.same_level_precursors(data, p)
                    for p in parent_candidates or possible_ions
                ]
            )

        # sometimes child peaks are heavier than parent
        return {k: v for k, v in children.items() if 0 < k < parent}

    def same_level_precursors(self, data, parent):
        result = {}
        adducts, tol = self.adduct_masses, self.tol
        for ion in data:
            target = parent - ion
            if any(
                d - tol < target + adduct < d + tol for d in data for adduct in adducts
            ):
                result[ion] = data[ion]
        return result
