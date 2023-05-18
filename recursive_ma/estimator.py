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


def upper_bound_by_MW(mw):
    return 0.06 * mw


def lower_bound_by_MW(mw):
    return log(max(mw * 0.06, 1.0), 2)


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
    def __init__(self, same_level=True, tol=0.01, exact_mass_digits=3, adduct_masses=COMMON_PRECURSORS):
        self.same_level = same_level
        self.tol = tol
        self.exact_mass_digits = exact_mass_digits
        self.isotope_weights = set(round(w, exact_mass_digits) for w in ISOTOPES.values())
        self.adduct_masses = adduct_masses

    def estimate_by_MW(self, mw):
        if round(mw, self.exact_mass_digits) in self.isotope_weights:
            print(f'HIT: {mw}')
            return interval([1.0, 1.0])
        return interval([lower_bound_by_MW(mw), upper_bound_by_MW(mw)])
    
    def estimate_MA(self, tree: dict[float, dict], mw: float, progress=False):
        children = tree.get(mw, None)
        mw_estimate = self.estimate_by_MW(mw)
        child_estimates = {}

        if not children:
            return mw_estimate
        else:
            for child in tqdm(children) if progress else children:
                complement = mw - child
                common_precursors = self.find_common_precursors(
                    children,
                    child,
                    complement,
                )

                ma_candidates = []
                for precursor in common_precursors | {0.0}:
                    chunks = [child - precursor, complement - precursor, precursor]
                    ma_candidates.append(
                        sum(
                            self.estimate_MA(
                                children,
                                chunk,
                            )
                            for chunk in chunks
                        )
                    )

                child_estimates[child] = min(ma_candidates, key=lambda x: x.midpoint)

            estimate = min(
                [mw_estimate, *child_estimates.values()], key=lambda x: x[0].sup
            )
            return estimate

    def find_common_precursors(self, data, parent1, parent2):
        precursors1 = self.precursors(data, parent1)
        precursors2 = self.precursors(data, parent2)
        return set(precursors1).intersection(precursors2)

    def precursors(self, data, parent):
        possible_ions = [parent + adduct for adduct in self.adduct_masses]
        parent_candidates = [
            d
            for d in data
            if any(d - self.tol < p < d + self.tol for p in possible_ions)
        ]
        children = unify_trees([
            {**{p - child: self.same_level_precursors(data, p - child) for child in data[p] or {}}, **(data[p] or {})}
            for p in parent_candidates
            ])
        if not children and self.same_level:
            children = unify_trees(
                [
                    self.same_level_precursors(data, p)
                    for p in parent_candidates or possible_ions
                ]
            )

        # sometimes child peaks are heavier than parent
        result = {k: v for k, v in children.items() if 0 < k < parent}
        return {**{p: data[p] for p in parent_candidates}, **result}

    def same_level_precursors(self, data, parent):
        result = {}
        tol = self.tol
        for ion in data:
            target = parent - ion
            if any(
                d - tol < target + adduct < d + tol
                for d in data
                for adduct in self.adduct_masses
            ):
                result[ion] = data[ion]
        return result
