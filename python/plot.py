import numpy as np
import matplotlib.pyplot as plt


training_performance = np.loadtxt('../cmake-build-debug/remote/output/4digits_1/training_performance.txt')
test_performance = np.loadtxt('../cmake-build-debug/remote/output/4digits_1/test_performance.txt')
training_performance_2 = np.loadtxt('../cmake-build-debug/remote/output/4digits_2/training_performance.txt')
test_performance_2 = np.loadtxt('../cmake-build-debug/remote/output/4digits_2/test_performance.txt')
title = "4 digits depth 0 vs depth 1"


# plt.plot(training_performance[:, 0], training_performance[:, 1], label='Training')
plt.plot(test_performance[:, 0], test_performance[:, 2], label='Test')
# plt.plot(training_performance_2[:, 0], training_performance_2[:, 1], label='Training depth 1')
plt.plot(test_performance_2[:, 0], test_performance_2[:, 2], label='Test depth 1')
plt.title('Accuracy ' + title)
plt.legend(loc='lower right')
plt.show()

# plt.plot(training_performance[:, 0], training_performance[:, 2], label='Training')
plt.plot(test_performance[:, 0], test_performance[:, 3], label='Test')
# plt.plot(training_performance_2[:, 0], training_performance_2[:, 2], label='Training depth 1')
plt.plot(test_performance_2[:, 0], test_performance_2[:, 3], label='Test depth 1')
plt.title('Error ' + title)
plt.legend(loc='upper right')
plt.show()

print('done')
