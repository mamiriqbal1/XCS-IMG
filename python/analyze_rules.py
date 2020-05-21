import numpy as np
import math

filter_size = 4  # 4x4

base_path = "../cmake-build-debug/output_2_6_2000_nokb_save_optimized/"

rules_file = open(base_path + "rule_with_codefragements.txt")


def parse_condition(line):
    values = []
    for comma_split in line.split(','):
        if comma_split.startswith('D0'):
            for space_split in comma_split.split(' '):
                if space_split != 'D0':
                    values.append(float(space_split))

    return values


def parse_action(line):
    value = int(line.split(" ")[1])
    return value


def parse_attributes(line):
    values = []
    tokens = line.split(' ')
    values.append((float(tokens[1])))
    values.append((float(tokens[3])))
    values.append((float(tokens[5])))
    values.append((float(tokens[7])))
    values.append((float(tokens[9])))
    values.append((float(tokens[11])))
    values.append((float(tokens[13])))
    values.append((float(tokens[15])))
    values.append((float(tokens[17])))
    values.append((float(tokens[19])))
    return values


def compare(window, filter_lower, filter_upper):
    assert (window.shape == filter_lower.shape)
    assert (window.shape == filter_upper.shape)
    if np.sum(window < filter_lower) > 0 or np.sum(window > filter_upper) > 0:
        return False
    return True


def does_match(img, filter):
    i_width = 28
    i_height = 28
    f_size = int(math.sqrt(filter.size/2))
    img = img.reshape(i_width, i_height)
    filter_lower = filter[::2].reshape(f_size, f_size)
    filter_upper = filter[1::2].reshape(f_size, f_size)

    for row in range(i_height - f_size + 1):
        for col in range(i_width - f_size + 1):
            if compare(img[row:row + f_size, col:col+f_size], filter_lower, filter_upper):
                return True
    return False


def count_matches_for_filter(good_filters, good_actions):
    mnist_3_8 = np.loadtxt('../data/mnist/mnist_train_3_8.txt')
    mnist_3_8 = mnist_3_8.round(2)
    for filter, action in zip(good_filters, good_actions):
        match_0 = 0
        match_1 = 0
        matched = 0
        for img in mnist_3_8:
            img_action = int(img[-1])
            img = img[:-1]
            if does_match(img, filter):
                matched += 1
                if img_action == 0:
                    match_0 += 1
                else:
                    match_1 += 1
        print('filter_action: '+ str(action) + '  matched: ' + str(matched) + '  action 0: ' + str(match_0) + '  action 1: ' + str(match_1))


filters = []
actions = []
attributes = []
line_number = 0
for line in rules_file:
    if line_number % 3 == 0:
        filters.append(parse_condition(line))
    elif line_number % 3 == 1:
        actions.append(parse_action(line))
    else:
        attributes.append(parse_attributes(line))
    line_number += 1

filters_np = np.array(filters)
actions_np = np.array(actions)
attributes_np = np.array(attributes)
fitness = attributes_np[:,3]
average_fitness = np.average(fitness)
good_index = fitness > average_fitness
good_filters = filters_np[good_index]
good_actions = actions_np[good_index]
good_attributes = attributes_np[good_index]
# np.savetxt(base_path + 'analyze/good_filters.txt', good_filters, delimiter=',')
# np.savetxt(base_path + 'analyze/good_actions.txt', good_actions, delimiter=',')
# np.savetxt(base_path + 'analyze/good_attributes.txt', good_attributes, delimiter=',')
# np.savetxt(base_path + 'analyze/all_filters.txt', filters_np, delimiter=',')
# np.savetxt(base_path + 'analyze/all_actions.txt', actions_np, delimiter=',')
# np.savetxt(base_path + 'analyze/all_attributes.txt', attributes_np, delimiter=',')
count_matches_for_filter(good_filters, good_actions)

print('done')
