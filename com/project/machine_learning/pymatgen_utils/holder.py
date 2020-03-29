from __future__ import division, print_function, unicode_literals

import math
import numpy as np
from numpy import array
# from warnings import warn
import sys

VERY_SMALL = 1E-12


# Common names for Holder means, by power:
# -1.0 harmonic mean
#  0.0 geometric mean
#  0.5 Hellinger mean
#  1.0 arithmetic mean
#  2.0 quadratic or Euclidean mean
#  3.0 cubic mean
#  4.0 quartic mean

def mean(values, power, weights=None, weights_norm=None):
    values = array(values, dtype=float)
    power = float(power)

    # Make sure weights match length and are normalized
    if weights is None:
        alpha = 1 / len(values)
    else:
        weights = array(weights, dtype=float)
        if len(values) != len(weights):
            # warn('Holder.mean returned zero when passed length mis-matched values and weights', UserWarning)
            sys.stderr.write('  Warning: Holder.mean returned zero when passed length mis-matched values and weights\n')
            return 0.0
        if weights_norm is not None:
            if weights_norm == "max" and max(weights) != 1.0:
                weights = weights / max(weights)
            elif weights_norm == "sum" and sum(weights) != 1.0:
                weights = weights / sum(weights)
            else:
                # warn('Holder.mean returned zero when passed unknown weights_norm method', UserWarning)
                sys.stderr.write('  Warning: Holder.mean returned zero when passed unknown weights_norm method\n')
                return 0.0
        alpha = 1 / sum(weights)

    if power == 0.0:  # geometric mean
        if any(abs(value) < VERY_SMALL for value in values):
            return 0.0
        elif any(value < 0 for value in values):
            # warn('Holder.mean returned zero when passed negative value with zero power', UserWarning)
            sys.stderr.write('  Warning: Holder.mean returned zero when passed negative value with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0

        if weights is None:
            return math.exp(alpha * sum(np.log(values)))
        else:
            return math.exp(alpha * sum(weights * np.log(values)))

    elif power == 1.0:  # arithmetic mean
        return np.average(values, weights=weights)

    if any(value < 0 for value in values):
        if power % 2 != 0.0:
            # warn('Holder.mean returned zero when passed negative value with non-even power', UserWarning)
            sys.stderr.write('  Warning: Holder.mean returned zero when passed negative value with non-even power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0

    if weights is None:
        return pow(alpha * sum(np.power(values, power)), 1 / power)
    else:
        return pow(alpha * sum(weights * np.power(values, power)), 1 / power)


def print_means(values, powers, weights=None, weights_norm=None):
    for power in powers:
        print(",{:.10f}".format(mean(values, power, weights, weights_norm)), end="")


def stdev(values, power, weights=None, weights_norm=None):
    values = array(values, dtype=float)
    power = float(power)

    # Check for single value in values
    if len(values) is 1:
        return 0.0

    # Make sure weights match length and are normalized
    if weights is None:
        beta = 1 / (len(values) - 1)
    else:
        weights = array(weights, dtype=float)
        if len(values) != len(weights):
            # warn('Holder.stdev returned zero when passed length mis-matched values and weights', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.stdev returned zero when passed length mis-matched values and weights\n')
            return 0.0
        if weights_norm is not None:
            if weights_norm == "max" and max(weights) != 1.0:
                weights = weights / max(weights)
            elif weights_norm == "sum" and sum(weights) != 1.0:
                weights = weights / sum(weights)
            else:
                # warn('Holder.stdev returned zero when passed unknown weights_norm method', UserWarning)
                sys.stderr.write('  Warning: Holder.stdev returned zero when passed unknown weights_norm method\n')
                return 0.0
        alpha = sum(weights)  # Note: Alpha is defined differently here than in mean function!
        beta = alpha / (alpha ** 2 - sum(np.power(weights, 2)))

    holder_mean = mean(values, power, weights=weights)

    if power == 0.0:  # geometric stdev (unbiased estimate)
        if any(value <= 0 for value in values):
            # warn('Holder.stdev returned zero when passed non-positive value with zero power', UserWarning)
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed non-positive value with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0
        if abs(holder_mean) < VERY_SMALL:
            # warn('Holder.stdev returned zero when passed values with near-zero mean with zero power', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.stdev returned zero when passed values with near-zero mean with zero power\n')
            sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
            return 0.0
        if weights is None:
            return math.exp(math.sqrt(beta * sum(np.power(np.log(values / holder_mean), 2))))
        else:
            return math.exp(math.sqrt(beta * sum(weights * np.power(np.log(values / holder_mean), 2))))

    holder_mean_centered_values = values - holder_mean

    if power == 1.0:  # arithmetic stdev (unbiased estimate)
        if weights is None:
            return math.sqrt(beta * sum(np.power(holder_mean_centered_values, 2)))
        else:
            return math.sqrt(beta * sum(weights * np.power(holder_mean_centered_values, 2)))

    # else:
    #   #warn('Holder.stdev returned zero when passed power other than 0 or 1', UserWarning)
    #   sys.stderr.write('  Warning: Holder.stdev returned zero when passed power other than 0 or 1\n')
    #   sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(values, weights, power))
    #   return 0.0

    # if sum(np.abs(holder_mean_centered_values)) < VERY_SMALL:
    #   return 0.0

    if power < 0.0:
        if any(abs(centered_value) < VERY_SMALL for centered_value in holder_mean_centered_values):
            # warn('Holder.stdev returned zero when passed near-zero value with negative power', UserWarning)
            sys.stderr.write('  Warning: Holder.stdev returned zero when passed near-zero value with negative power\n')
            sys.stderr.write('    values = {:s}  power = {:f}\n'.format(holder_mean_centered_values, power))
            return 0.0

    if any(centered_value < 0 for centered_value in holder_mean_centered_values):
        if power % 1 != 0.0:
            # warn('Holder.stdev returned zero when passed non-positive value with non-integer power', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.stdev returned zero when passed negative value with non-integer power\n')
            sys.stderr.write('    values = {:s}  power = {:f}\n'.format(holder_mean_centered_values, power))
            return 0.0

    if weights is None:
        # sys.stderr.write('    values = {:s}  power = {:f}\n'.format(holder_mean_centered_values, power))
        return pow(beta * sum(np.power(holder_mean_centered_values, 2 * power)), 1 / 2 / power)
    else:
        # sys.stderr.write('    values = {:s}  weights = {:s}  power = {:f}\n'.format(holder_mean_centered_values, weights, power))
        return pow(beta * sum(weights * np.power(holder_mean_centered_values, 2 * power)), 1 / 2 / power)


def print_stdevs(values, powers, weights=None, weights_norm=None):
    for power in powers:
        print(",{:.10f}".format(stdev(values, power, weights, weights_norm)), end="")


def cfvar(values, power, weights=None, weights_norm=None):
    values = array(values, dtype=float)
    power = float(power)
    length = len(values)

    # Make sure weights match length and are normalized
    if weights is not None:
        weights = array(weights, dtype=float)
        if length != len(weights):
            # warn('Holder.cfvar returned zero when passed length mis-matched values and weights', UserWarning)
            sys.stderr.write(
                '  Warning: Holder.cfvar returned zero when passed length mis-matched values and weights\n')
            return 0.0
        if weights_norm is not None:
            if weights_norm == "max" and max(weights) != 1.0:
                weights = weights / max(weights)
            elif weights_norm == "sum" and sum(weights) != 1.0:
                weights = weights / sum(weights)
            else:
                # warn('Holder.cfvar returned zero when passed unknown weights_norm method', UserWarning)
                sys.stderr.write('  Warning: Holder.cfvar returned zero when passed unknown weights_norm method\n')
                return 0.0

    holder_mean = mean(values, power, weights=weights)
    holder_stdev = stdev(values, power, weights=weights)

    if holder_stdev == 0:
        return 0.0  # Regardless of mean
    elif holder_mean == 0:
        # warn('Holder.cfvar returned zero when passed values with zero mean', UserWarning)
        sys.stderr.write('  Warning: Holder.cfvar returned zero when passed values with zero mean\n')
        sys.stderr.write('    values = {:s}  power = {:f}\n'.format(values, power))
        return 0.0
    else:
        return holder_stdev / holder_mean


def print_cfvars(values, powers, weights=None, weights_norm=None):
    for power in powers:
        print(",{:.10f}".format(cfvar(values, power, weights, weights_norm)), end="")
