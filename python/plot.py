import numpy as np
import matplotlib.pyplot as plt


training_performance = np.loadtxt('../cmake-build-debug/remote/output/2digits_1/training_performance.txt')
test_performance = np.loadtxt('../cmake-build-debug/remote/output/2digits_1/test_performance.csv.txt')


plt.plot(training_performance[:, 0], training_performance[:, 1])
plt.plot(test_performance[:, 0], test_performance[:, 2])
plt.title('Training Accuracy digits 3, 8')
plt.legend(['Training Accuracy digits 3, 8'])
plt.show()

plt.plot(training_performance[:, 0], training_performance[:, 2])
plt.plot(test_performance[:, 0], test_performance[:, 3])
plt.title('Training Error digits 3, 8')
plt.legend(['Training Error digits 3, 8'])
plt.show()


print('done')