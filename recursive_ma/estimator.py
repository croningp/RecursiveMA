from math import log
from interval import interval
from functools import reduce


def upper_bound_by_MW(mw):
    return 0.05 * mw


def lower_bound_by_MW(mw):
    return log(max(mw / 20.0, 1.0), 2)


def estimate_by_MW(mw):
    return interval([lower_bound_by_MW(mw), upper_bound_by_MW(mw)])


def estimate_MA(data, mw, same_level=True, decimals=1):
    children = data.get(mw, None)
    mw_estimate = estimate_by_MW(mw)

    if not children:
        return mw_estimate
    else:
        child_estimates = {
            child: estimate_MA(
                children, child, same_level=same_level, decimals=decimals
            )
            for child in children
        }
        for child in children:
            # TODO: Use same number of decimal places as original data
            complement = round(mw - child, decimals)
            child_estimates[child] += estimate_MA(
                children,
                complement,
                same_level=same_level,
                decimals=decimals,
            )
            child_estimates[child] += 1  # one composition step
            child_estimates[child] -= common_MA(
                children,
                child,
                complement,
                same_level=same_level,
                decimals=decimals,
            )
        # intersection corrected estimates from children and self
        estimate = reduce(lambda x, y: x & y, [mw_estimate, *child_estimates.values()])
        return estimate


def same_level_precursors(data, parent, decimals):
    result = {}
    for ion in data:
        # TODO: Use same number of decimal places as original data
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


def common_MA(data, parent1, parent2, same_level, decimals):
    precursors1 = precursors(data, parent1, same_level=same_level, decimals=decimals)
    precursors2 = precursors(data, parent2, same_level=same_level, decimals=decimals)
    common_precursors = set(precursors1).intersection(precursors2)
    ion_MAs = [
        estimate_MA(precursors1, precursor, same_level=same_level, decimals=decimals)
        & estimate_MA(precursors2, precursor, same_level=same_level, decimals=decimals)
        for precursor in common_precursors
    ]
    return (
        max(ion_MAs, key=lambda ma_interval: ma_interval.midpoint)
        if ion_MAs
        else interval([0.0, 0.0])
    )
