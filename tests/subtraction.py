"""Tools for computing what cmuts-sub should write from its input files.

Each rule is applied in the type its field is stored as, matching cmuts-sub:
rates in float32, counts as whole unsigneds. A sum or a difference is therefore
exact, while a ratio and its error round enough that a caller should allow a
tolerance.
"""

from __future__ import annotations

import numpy as np

from outputs import COVERAGE, ERROR, LENGTHS, COUNTED, REACTIVITY, REJECTED, UNMAPPED
from outputs import field_of

# Every rule takes the values of the field it forms and the reactivities of the
# same inputs, one array apiece. Only the error of a ratio reads a field other
# than its own.


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
    return _over_control(values[0] - values[1], values[2])


def _ratio_error(values, rates):
    ratio = _over_control(rates[0] - rates[1], rates[2])
    spread = values[0] ** 2 + values[1] ** 2 + (ratio * values[2]) ** 2

    return _over_control(np.sqrt(spread), rates[2])


# How each dataset of the output is formed from the inputs, without a denatured
# control and with one.
RULES = {
    COVERAGE: (_add, _add),
    REACTIVITY: (_subtract, _ratio),
    ERROR: (_quadrature, _ratio_error),
    LENGTHS: (_add, _add),
    COUNTED: (_add, _add),
    REJECTED: (_add, _add),
    UNMAPPED: (_add, _add),
}

CONTROLLED = 1
UNCONTROLLED = 0


def expected(name, *inputs):
    """Computes the values a single field of the output should hold.

    Takes whole input files, since the error of a ratio reads the reactivities
    as well as the errors.
    """
    rule = RULES[name][CONTROLLED if len(inputs) > 2 else UNCONTROLLED]
    values = [np.asarray(field_of(path, name)) for path in inputs]
    rates = [np.asarray(field_of(path, REACTIVITY)) for path in inputs]

    return rule(values, rates)
