# Abhishek/Michael MW => MS formula (upper bound)
def estimate_by_MW(mw):
    return 0.05 * mw + 2.5


def estimate_MA(data, mw):
    children = data[mw]
    mw_estimate = estimate_by_MW(mw)

    if children == {}:
        # does not fragment: assume single atom
        # return 1.0
        return mw_estimate
    elif children is None:
        # haven't tried to fragment this ion: use MW estimate
        return mw_estimate
    else:
        child_estimates = {child: estimate_MA(children, child) for child in children}
        for child in children:
            # attempt to find complementary ion
            complement = min(children, key=lambda c: abs(c + child - mw))
            if abs(child + complement - mw) < 2.0:
                child_estimates[child] += estimate_MA(children, complement)
                child_estimates[child] -= common_MA(children, child, complement)
            else:
                child_estimates[child] += estimate_by_MW(mw - child)
        # min corrected estimates from children and self
        return min([mw_estimate, *child_estimates.values()])


def same_level_precursors(data, parent):
    result = []
    for ion in data:
        if parent - ion in data:
            result.append(ion)
    return set(result)


def common_MA(data, parent1, parent2):
    if not data[parent1] or not data[parent2]:
        # estimate using precursors from the same level
        common_ions = same_level_precursors(data, parent1).intersection(
            same_level_precursors(data, parent2)
        )
        ion_MAs = [estimate_MA(data, ion) for ion in common_ions]
        return max(ion_MAs) if ion_MAs else 0.0

    common_ions = [
        (ion1, ion2)
        for ion1 in data[parent1]
        for ion2 in data[parent2]
        if abs(ion1 - ion2) < 0.1
    ]

    ion_MAs = [
        min(estimate_MA(data[parent1], ion1), estimate_MA(data[parent2], ion2))
        for ion1, ion2 in common_ions
    ]

    return max(ion_MAs) if ion_MAs else 0.0
