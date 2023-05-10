from math import log
from interval import interval
from functools import reduce

from tqdm.auto import tqdm


def upper_bound_by_MW(mw):
    return 0.05 * mw


def lower_bound_by_MW(mw):
    return log(max(mw / 20.0, 1.0), 2)


def estimate_by_MW(mw):
    return interval([lower_bound_by_MW(mw), upper_bound_by_MW(mw)])


class MAEstimator:
    def __init__(self, same_level=True, decimals=1, adduct_mass=1.00):
        self.same_level = same_level
        self.decimals = decimals
        # default adduct is H+ (monoisotopic mass = 1.007825)
        self.adduct_mass = adduct_mass

    def estimate_MA(self, tree: dict[float, dict], mw: float, progress=False):
        children = tree.get(mw, None)
        mw_estimate = estimate_by_MW(mw)
        child_estimates = {}

        if not children:
            return mw_estimate
        else:
            for child in tqdm(children) if progress else children:
                complement = round(mw - child + self.adduct_mass, self.decimals)
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
                                round(chunk, self.decimals),
                            )
                            for chunk in chunks
                        )
                    )

                child_estimates[child] = reduce(lambda x, y: x & y, ma_candidates)

            # intersection corrected estimates from children and self
            estimate = reduce(
                lambda x, y: x & y, [mw_estimate, *child_estimates.values()]
            )
            return estimate

    def find_common_precursors(self, data, parent1, parent2):
        precursors1 = self.precursors(data, parent1)
        precursors2 = self.precursors(data, parent2)
        return set(precursors1).intersection(precursors2)

    def precursors(self, data, parent):
        children = data.get(parent, None)
        if not children:
            result = self.same_level_precursors(data, parent) if self.same_level else {}
        else:
            # sometimes child peaks are heavier than parent
            result = {k: v for k, v in children.items() if k < parent}
        return {parent: children, **result}

    def same_level_precursors(self, data, parent):
        result = {}
        for ion in data:
            if round(parent - ion + self.adduct_mass, self.decimals) in data:
                result[ion] = None
        return result
