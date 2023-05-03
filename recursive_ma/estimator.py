# Abhishek/Michael MW => MS formula (upper bound)
def estimate_by_MW(mw):
    return 0.05 * mw + 2.5


def estimate_MA(data, mw):
    children = data.get(mw, None)
    mw_estimate = estimate_by_MW(mw)

    if not children or mw not in children:
        return mw_estimate
    else:
        child_estimates = {child: estimate_MA(children, child) for child in children}
        for child in children:
            # TODO: Use same number of decimal places as original data
            complement = round(mw - child, 1)
            child_estimates[child] += estimate_MA(children, complement)
            child_estimates[child] += 1 # one composition step
            child_estimates[child] -= common_MA(children, child, complement)
        # min corrected estimates from children and self
        improved_estimates = [estimate for estimate in child_estimates.values() if estimate < mw_estimate] + [mw_estimate]
        return sorted([mw_estimate, *improved_estimates])[len(improved_estimates) // 2]


def same_level_precursors(data, parent):
    result = {}
    for ion in data:
        # TODO: Use same number of decimal places as original data
        if round(parent - ion, 1) in data:
            result[ion] = None
    return result


def precursors(data, parent):
    if not data.get(parent, None):
        return same_level_precursors(data, parent)
    else:
        return data[parent]


def common_MA(data, parent1, parent2):
    precursors1 = precursors(data, parent1)
    precursors2 = precursors(data, parent2)
    common_precursors = set(precursors1).intersection(precursors2)
    ion_MAs = [
        min(estimate_MA(precursors1, precursor), estimate_MA(precursors2, precursor))
        for precursor in common_precursors
    ]
    return max(ion_MAs) if ion_MAs else 0.0
