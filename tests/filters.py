"""The criteria a run is given, and how to name one in a test report."""

# Every criterion at once. The length bounds read the stored sequence, which a
# CRAM reconstructs on decode, so this is the set to try where the format may
# decide the answer.
COMPOUND = {"min_mapq": 30, "strand": "reverse", "min_length": 100, "max_length": 500}

# Combinations chosen so that each criterion is exercised on its own, at a
# boundary, and alongside the others.
FILTERS = [
    {},
    {"min_mapq": 1},
    {"min_mapq": 30},
    {"min_mapq": 60},
    {"strand": "forward"},
    {"strand": "reverse"},
    {"min_length": 200},
    {"max_length": 300},
    {"min_length": 150, "max_length": 400},
    {"min_length": 9000},                      # accepts no read
    COMPOUND,
]

# Criteria are given in full at every call, so no default is relied on.
UNFILTERED = {"min_mapq": 0, "strand": "forward,reverse", "min_length": 0, "max_length": 0}


def criteria(filters: dict) -> dict:
    """Returns a full set of criteria, with the given filters applied over the
    unfiltered defaults."""
    return {**UNFILTERED, **filters}


def describe_filters(filters: dict) -> str:
    """Returns the name a set of filters is reported under."""
    return ",".join(f"{k}={v}" for k, v in filters.items()) or "unfiltered"
