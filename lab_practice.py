from sympy import Symbol, solve
import numpy as np
from math import log2, floor

def generate_sum(coll):
    sum = ''
    for x in coll:
        if (str(x)[0] == '-'):
            pass
        else:
            sum += '+'
        sum += f'{x}'

    return sum

def create_transform_matrix(q_filters, p_filters, n):
    A = []

    for i in range(n // 2):
        k = (2 * i) % n
        q_vector = []
        p_vector = []
        null_vector = [0 for x in range(k)]
        offset = k + len(q_filters) - n

        if offset <= 0:
            q_vector = null_vector + q_filters + [0 for x in range(abs(offset))]
            p_vector = null_vector + p_filters + [0 for x in range(abs(offset))]
        else:
            q_vector = q_filters[-offset:] + [0 for x in range(k - offset)] + q_filters[:-offset]
            p_vector = p_filters[-offset:] + [0 for x in range(k - offset)] + p_filters[:-offset]

        A.append(q_vector)
        A.append(p_vector)

    return A

def transform_signal(A, s):
    new_signal = np.dot(A, s)
    hi = [new_signal[2 * i] for i in range(len(new_signal) // 2)]
    lo = [new_signal[2 * i + 1] for i in range(len(new_signal) // 2)]

    return hi + lo

def mid_level_signal_transform(q, p, s):
    A = create_transform_matrix(q, p, len(s))
    result = [[transform_signal(A, s)]]

    for i in range(floor(log2(len(s))) // 2):
        A = create_transform_matrix(q, p, len(result[i][0]) // 2)
        level_transform = []
        
        for j in range(len(result[i])):
            x = result[i][j]
            level_transform.append(transform_signal(A, x[:len(x) // 2]))
            level_transform.append(transform_signal(A, x[len(x) // 2:]))

        result.append(level_transform)

    return result[-1]

def formulate_hi_fr_filters(n):
    knowledge_coeff = [Symbol(f'K{x}') for x in range(n * (n // 2))]
    hi_fr_filters = [Symbol(f'q{x}') for x in range(n)]

    eq1 = [x ** 2 for x in hi_fr_filters]
    eq1.append(-1)

    eq1 = eq1
    null_moments_eqs = [[hi_fr_filters[x] * (-(n // 2) + 1 + x) ** b for x in range(n)] for b in range(n // 2)]    
    class_restrictions_eqs = [[hi_fr_filters[x] * knowledge_coeff[i*n + x] for x in range(n)] for i in range(n // 2 - 1)]

    non_linear_eq_system_1 = [generate_sum(eq1)]
    for x in null_moments_eqs:
        non_linear_eq_system_1.append(generate_sum(x))
    for x in class_restrictions_eqs:
        non_linear_eq_system_1.append(generate_sum(x))

        
    sol = solve(non_linear_eq_system_1, hi_fr_filters)

    return sol

def formulate_lo_fr_filters(q_filters):
    return [((-1) ** (i + 1)) * q_filters[len(q_filters) - i - 1] for i in range(len(q_filters))]

# q_filters = formulate_hi_fr_filters(4)
# print(q_filters)
# p_filters = formulate_lo_fr_filters(q_filters)

# A = np.array(create_transform_matrix(q_filters, p_filters, len(s_input)))
# s_input = np.array(s_input)

# transform_data = transform_signal_data(A, s)
s_input = [3, 4, 5, 8, 2, 1, 2, 4, 2, 5, 6, 7, 5, 3, 4, 5, 3, 4, 5, 8, 2, 1, 2, 4, 2, 5, 6, 7, 5, 3, 4, 5, 3, 4, 5, 8, 2, 1, 2, 4, 2, 5, 6, 7, 5, 3, 4, 5, 3, 4, 5, 8, 2, 1, 2, 4, 2, 5, 6, 7, 5, 3, 4, 5]
A = create_transform_matrix([1, 4, 5, 2, 3, 6], [5, 4, 6, 7, 8, 9], len(s_input))
# s = transform_signal_data(A, s_input, len(s_input))
# print(np.shape(s))
s = mid_level_signal_transform([1, 4, 5, 2, 3, 6], [5, 4, 6, 7, 8, 9], s_input)
print(len(s))

# print(q_filters)
# print(p_filters)
