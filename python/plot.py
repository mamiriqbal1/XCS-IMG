import numpy as np
import matplotlib.pyplot as plt


training_performance = np.loadtxt('../cmake-build-debug/remote/output/4digits_2/training_performance.txt')
test_performance = np.loadtxt('../cmake-build-debug/remote/output/4digits_2/test_performance.txt')
title = "4 digits, depth=1"


plt.plot(training_performance[:, 0], training_performance[:, 1], label='Training')
plt.plot(test_performance[:, 0], test_performance[:, 2], label='Test')
plt.title('Training Accuracy ' + title)
plt.legend(loc='lower right')
plt.show()

plt.plot(training_performance[:, 0], training_performance[:, 2], label='Training')
plt.plot(test_performance[:, 0], test_performance[:, 3], label='Test')
plt.title('Training Error ' + title)
plt.legend(loc='upper right')
plt.show()

print('done')
