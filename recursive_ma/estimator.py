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


def estimate_MA(data, mw, same_level=True, decimals=1, progress=False):
    children = data.get(mw, None)
    mw_estimate = estimate_by_MW(mw)
    child_estimates = {}

    if not children:
        return mw_estimate
    else:
        for child in tqdm(children) if progress else children:
            complement = round(mw - child, decimals)
            precursors = common_precursors(
                children,
                child,
                complement,
                same_level=same_level,
                decimals=decimals,
            )

            ma_candidates = []
            for precursor in precursors | {0.0}:
                chunks = [child - precursor, complement - precursor, precursor]
                ma_candidates.append(
                    sum(
                        estimate_MA(
                            children,
                            round(chunk, decimals),
                            same_level=same_level,
                            decimals=decimals,
                        )
                        for chunk in chunks
                    )
                )

            child_estimates[child] = reduce(lambda x, y: x & y, ma_candidates)

        # intersection corrected estimates from children and self
        estimate = reduce(lambda x, y: x & y, [mw_estimate, *child_estimates.values()])
        return estimate


def same_level_precursors(data, parent, decimals):
    result = {}
    for ion in data:
        if round(parent - ion, decimals) in data:
            result[ion] = None
    return result


def precursors(data, parent, same_level, decimals):
    children = data.get(parent, None)
    if not children:
        result = same_level_precursors(data, parent, decimals) if same_level else {}
    else:
        # sometimes child peaks are heavier than parent
        result = {k: v for k, v in children.items() if k < parent}
    return {parent: children, **result}


def common_precursors(data, parent1, parent2, same_level, decimals):
    precursors1 = precursors(data, parent1, same_level=same_level, decimals=decimals)
    precursors2 = precursors(data, parent2, same_level=same_level, decimals=decimals)
    return set(precursors1).intersection(precursors2)
