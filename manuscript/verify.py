import json
from fractions import Fraction
from pathlib import Path
from math import comb
from math import isqrt


def mobius_trial(number):
    remaining = number
    sign = 1
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            remaining //= divisor
            sign = -sign
            if remaining % divisor == 0:
                return 0
        divisor += 1
    return -sign if remaining > 1 else sign


def mobius_sieve(limit):
    values = [1] * (limit + 1)
    primes = bytearray(b'\x01') * (limit + 1)
    for candidate in range(2, limit + 1):
        if not primes[candidate]:
            continue
        for multiple in range(candidate, limit + 1, candidate):
            values[multiple] *= -1
            if multiple != candidate:
                primes[multiple] = 0
        square = candidate * candidate
        for multiple in range(square, limit + 1, square):
            values[multiple] = 0
    values[0] = 0
    return values


def is_squarefree(number):
    return all(number % (divisor * divisor) for divisor in range(2, isqrt(number) + 1))


def sequence_value(number, mu):
    return 0 if number % 2 else -mu


def collatz_step(number):
    return 3 * number + 1 if number % 2 else number // 2


def require(condition, label, number):
    if not condition:
        raise ValueError(f'{label}: n={number}')


def divisors(number):
    return [candidate for candidate in range(1, number + 1) if number % candidate == 0]


def prime_exponents(number):
    factors = []
    candidate = 2
    while candidate * candidate <= number:
        exponent = 0
        while number % candidate == 0:
            exponent += 1
            number //= candidate
        if exponent:
            factors.append((candidate, exponent))
        candidate += 1
    if number > 1:
        factors.append((number, 1))
    return factors


def root_local(exponent):
    value = Fraction(1)
    for index in range(exponent):
        value *= -(Fraction(1, 2) - index) / (index + 1)
    return value


def outward_interval(lower, upper):
    scale = 10 ** 12
    lower_integer = lower.numerator * scale // lower.denominator
    upper_integer = -((-upper.numerator * scale) // upper.denominator)
    return {'lower_times_10_pow_12': lower_integer,
            'upper_times_10_pow_12': upper_integer}


def sibling_checks(mu, limit):
    convolution_limit = 2048
    root_limit = 512
    root = [Fraction(0)] * (root_limit + 1)
    root[1] = Fraction(1)
    for number in range(1, convolution_limit + 1):
        divs = divisors(number)
        companion = (number % 2) * mu[number]
        power_two = int(number & (number - 1) == 0)
        inverse = sum((divisor % 2) * mu[divisor] * ((number // divisor) % 2)
                      for divisor in divs)
        transformed = sum(mu[number // divisor] for divisor in divs
                          if divisor & (divisor - 1) == 0)
        inverse_transformed = sum((divisor % 2) * mu[divisor] for divisor in divs)
        require(inverse == int(number == 1), 'dirichlet_inverse', number)
        require(transformed == companion, 'mobius_transform', number)
        require(inverse_transformed == power_two, 'inverse_mobius_transform', number)
        require(-mobius_trial(2 * number) == companion, 'companion_halving', number)
        discrepancy = companion != mu[number]
        require(discrepancy == (number % 2 == 0 and mu[number] != 0),
                'empirical_deviation', number)
        transformed_mu = mu[number] - (mu[number // 2] if number % 2 == 0 else 0)
        require((abs(transformed_mu) == 2) == (number % 2 == 0 and mu[number] != 0),
                'A092673_support', number)
        if number <= root_limit:
            if number > 1:
                root[number] = (Fraction(companion) - sum(
                    (root[divisor] * root[number // divisor]
                     for divisor in divs if 1 < divisor < number), Fraction(0))) / 2
            expected = Fraction(1) if number % 2 else Fraction(0)
            denominator = 1
            for prime, exponent in prime_exponents(number):
                expected *= root_local(exponent)
                denominator *= Fraction(comb(2 * exponent, exponent), 4 ** exponent).denominator
            require(root[number] == expected, 'root_euler_coefficients', number)
            require(sum(root[divisor] * root[number // divisor] for divisor in divs) == companion,
                    'rational_convolution_square', number)
            if number % 2:
                require(root[number].denominator == denominator, 'odd_root_denominator', number)
    require(mu[27] == 0 and (Fraction(-1) + (-1) ** 3) / 2 == -1,
            'prime_power_correction', 27)
    generalized_steps = 0
    for multiplier, shift in [(1, 1), (3, 1), (5, 3), (7, 5)]:
        for start in range(1, 101):
            current = start
            accumulated = 0
            initial = sequence_value(start, mobius_trial(start))
            for index in range(12):
                following = multiplier * current + shift if current % 2 else current // 2
                accumulated += mobius_trial(following)
                require(accumulated == initial - sequence_value(following, mobius_trial(following)),
                        'generalized_trajectory', start)
                current = following
                generalized_steps += 1
    require(mobius_trial(5) != sequence_value(3, -1) - sequence_value(5, -1),
            'accelerated_map_counterexample', 3)
    rational_constants = {}
    cutoff = 256
    for exponent in [2, 3, 4]:
        companion_sum = sum((Fraction((number % 2) * mu[number], number ** exponent)
                             for number in range(1, cutoff + 1)), Fraction(0))
        odd_sum = sum((Fraction(1, number ** exponent) for number in range(1, cutoff + 1, 2)), Fraction(0))
        zeta_sum = sum((Fraction(1, number ** exponent) for number in range(1, cutoff + 1)), Fraction(0))
        tail = Fraction(1, (exponent - 1) * cutoff ** (exponent - 1))
        lower, upper = companion_sum - tail, companion_sum + tail
        inverse_lower, inverse_upper = 1 / (odd_sum + tail), 1 / odd_sum
        zeta_factor = 1 - Fraction(1, 2 ** exponent)
        product_lower = 1 / (zeta_factor * (zeta_sum + tail))
        product_upper = 1 / (zeta_factor * zeta_sum)
        require(max(lower, inverse_lower, product_lower) <= min(upper, inverse_upper, product_upper),
                'dirichlet_constant_intervals', exponent)
        if exponent in [3, 4]:
            require(upper < odd_sum, 'nonreciprocal_constant_refuted', exponent)
        rational_constants[str(exponent)] = {
            'B_s': outward_interval(lower, upper),
            'reciprocal_odd_power_sum': outward_interval(inverse_lower, inverse_upper),
            'odd_power_sum': outward_interval(odd_sum, odd_sum + tail),
        }
    support = [number for number in range(1, limit + 1) if number % 2 == 0 and mu[number]]
    odd_support = [number for number in range(1, limit + 1) if number % 2 and mu[number]]
    counts = {}
    for prime in [2, 3, 5, 7]:
        count = sum(mu[number] != 0 and number % prime == 0 for number in range(1, limit + 1))
        coprime_count = sum(mu[number] != 0 and number % prime != 0 for number in range(1, limit // prime + 1))
        require(count == coprime_count, 'prime_support_bijection', prime)
        identity_count = sum(mu[divisor] * (limit // (divisor * divisor)
                            - limit // (prime * divisor * divisor))
                            for divisor in range(1, isqrt(limit) + 1) if divisor % prime)
        require(identity_count == sum(mu[number] != 0 and number % prime != 0
                                     for number in range(1, limit + 1)),
                'squarefree_counting_formula', prime)
        counts[str(prime)] = count
    records = json.loads(Path(__file__).with_name('source-values.json').read_text())
    checked_values = {}
    for sequence, record in records.items():
        for index, expected in enumerate(record['values'], record['offset']):
            if sequence == 'A292273': actual = sequence_value(index, mu[index])
            elif sequence == 'A087003': actual = (index % 2) * mu[index]
            elif sequence == 'A039956': actual = support[index - 1]
            elif sequence == 'A056911': actual = odd_support[index - 1]
            elif sequence == 'A264387': actual = (support[index - 1] - 2) // 4
            elif sequence == 'A099990': actual = mu[2 * index + 1]
            elif sequence == 'A099991': actual = mu[2 * index]
            elif sequence == 'A209229': actual = int(index > 0 and index & (index - 1) == 0)
            elif sequence == 'A318657': actual = root[index].numerator
            elif sequence == 'A318658': actual = root[index].denominator
            elif sequence == 'A046644':
                actual = 1
                for prime, exponent in prime_exponents(index):
                    actual *= Fraction(comb(2 * exponent, exponent), 4 ** exponent).denominator
            elif sequence == 'A092673': actual = mu[index] - (mu[index // 2] if index % 2 == 0 else 0)
            elif sequence == 'A091069': actual = int(index % 2 != 0 and mu[index] != 0) * (1 if index % 8 in [1, 7] else -1)
            else: raise ValueError(sequence)
            require(actual == expected, 'source_' + sequence, index)
        checked_values[sequence] = len(record['values'])
    return {
        'divisor_identity_n_max': convolution_limit,
        'rational_square_root_n_max': root_limit,
        'prime_power_counterexample': {'n': 27, 'source_expression': -1, 'correct_value': 0},
        'accelerated_map_counterexample': {'segment': [3, 5], 'mobius_sum': -1, 'endpoint_difference': 0},
        'generalized_map_prefix_checks': generalized_steps,
        'generalized_map_odd_branches': ['n+1', '3n+1', '5n+3', '7n+5'],
        'constant_interval_cutoff': cutoff,
        'constant_rational_intervals': rational_constants,
        'squarefree_multiple_counts_at_65537': counts,
        'displayed_source_values_checked': checked_values,
    }


def main():
    limit = 65537
    crosscheck_limit = 10000
    trajectory_limit = 1000
    trajectory_step_cap = 10000
    mu = mobius_sieve(limit)
    source_values = [
        0, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0,
        0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0,
        0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0,
        0, 1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0,
        0, -1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0,
    ]
    real_cycle = (1, 0, -1, 0)
    residue_counts = [0, 0, 0, 0]
    support_count = 0
    for number in range(1, limit + 1):
        value = sequence_value(number, mu[number])
        residue_counts[number % 4] += 1
        require(value == mu[number] * real_cycle[number % 4], 'real_part', number)
        require(value == (number % 2 - 1) * mu[number], 'parity', number)
        require(value == (number % 2) * mu[number] - mu[number], 'companion', number)
        if number % 2 == 0:
            half = number // 2
            require(value == (mu[half] if half % 2 else 0), 'halving', number)
        require((value != 0) == (number % 2 == 0 and is_squarefree(number)), 'support', number)
        support_count += value != 0
    for number in range(1, crosscheck_limit + 1):
        require(mu[number] == mobius_trial(number), 'mobius_crosscheck', number)
    for number, expected in enumerate(source_values, 1):
        require(sequence_value(number, mu[number]) == expected, 'source_values', number)
    max_steps = 0
    max_value = 1
    total_steps = 0
    for start in range(1, trajectory_limit + 1):
        current = start
        accumulated = 0
        steps = 0
        initial = sequence_value(start, mobius_trial(start))
        while current != 1:
            require(steps < trajectory_step_cap, 'trajectory_bound', start)
            following = collatz_step(current)
            following_mu = mobius_trial(following)
            current_mu = mobius_trial(current)
            current_b = (current % 2) * current_mu
            following_b = (following % 2) * following_mu
            require(current_b - following_b == current_mu, 'coboundary', current)
            accumulated += following_mu
            current = following
            steps += 1
            max_value = max(max_value, current)
            require(accumulated == initial - sequence_value(current, following_mu), 'finite_trajectory', start)
        require(accumulated == initial, 'terminal_sum', start)
        require(accumulated + mobius_trial(start) == (start % 2) * mobius_trial(start), 'inclusive_sum', start)
        total_steps += steps
        max_steps = max(max_steps, steps)
    siblings = sibling_checks(mu, limit)
    print(json.dumps({
        'status': 'pass',
        'scope': 'Finite exact integer checks; all-index proofs are in note.tex.',
        'identity_n_max': limit,
        'mobius_trial_crosscheck_n_max': crosscheck_limit,
        'source': 'OEIS A292273 revision 24, displayed initial values, offset 1',
        'source_values_checked': len(source_values),
        'residue_counts_0_1_2_3': residue_counts,
        'nonzero_positions_checked': support_count,
        'trajectory_start_max': trajectory_limit,
        'trajectory_step_cap': trajectory_step_cap,
        'trajectory_total_steps': total_steps,
        'trajectory_max_steps': max_steps,
        'trajectory_max_value': max_value,
        'sibling_checks': siblings,
        'checks': ['real_part', 'parity', 'companion', 'halving', 'support',
                   'mobius_crosscheck', 'source_values', 'coboundary',
                   'finite_trajectory', 'terminal_sum', 'inclusive_sum'],
    }, ensure_ascii=False, indent=2) + '\n', end='')


if __name__ == '__main__':
    main()
