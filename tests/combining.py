"""Tools for computing what a combination should write from its input files.

Each rule is applied in the type its field is stored as, matching the programs:
rates in float32, counts as whole unsigneds. A sum or a difference is therefore
exact, while a ratio and its error round enough that a caller should allow a
tolerance.
"""

from __future__ import annotations

import numpy as np

from outputs import COVERAGE, ERROR, LENGTHS, COUNTED, REACTIVITY, REJECTED, UNMAPPED
from outputs import field_of

# Every rule takes the values of the field being formed and the reactivities of
# the same inputs, one array apiece. Only the error of a ratio uses a field
# other than the one being formed.


def _add(values, rates):
    total = values[0]

    for value in values[1:]:
        total = total + value

    return total


def _subtract(values, rates):
    return values[0] - values[1]


def _quadrature(values, rates):
    return np.sqrt(values[0] * values[0] + values[1] * values[1])


def _over_control(values, denatured):
    """Divides by the denatured control. The result is NaN wherever the
    control's reactivity is not above zero."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denatured > 0, values / denatured, np.float32(np.nan))


def _ratio(values, rates):
    return _over_control(values[0], rates[1])


def _ratio_error(values, rates):
    ratio = _over_control(rates[0], rates[1])
    spread = values[0] ** 2 + (ratio * values[1]) ** 2

    return _over_control(np.sqrt(spread), rates[1])


# How each dataset of the output is formed from the inputs, for each program.
# Both sum every count, so only the rate and its error differ.
SUB_RULES = {
    COVERAGE: _add,
    REACTIVITY: _subtract,
    ERROR: _quadrature,
    LENGTHS: _add,
    COUNTED: _add,
    REJECTED: _add,
    UNMAPPED: _add,
}

DIV_RULES = {
    COVERAGE: _add,
    REACTIVITY: _ratio,
    ERROR: _ratio_error,
    LENGTHS: _add,
    COUNTED: _add,
    REJECTED: _add,
    UNMAPPED: _add,
}


def expected(rules, name, *inputs):
    """Computes the values a single field of the output should hold.

    Takes whole input files, since the error of a ratio reads the reactivities
    as well as the errors.
    """
    values = [np.asarray(field_of(path, name)) for path in inputs]
    rates = [np.asarray(field_of(path, REACTIVITY)) for path in inputs]

    return rules[name](values, rates)
