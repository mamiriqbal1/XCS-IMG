import numpy as np
import matplotlib.pyplot as plt


training_performance = np.loadtxt('../cmake-build-debug/output_4_3/output_training_.txt')
training_performance_nokb = np.loadtxt('../cmake-build-debug/output_4_3_nokb/output_training_.txt')
training_performance_nokb_2 = np.loadtxt('../cmake-build-debug/output_2_3_nokb/output_training_.txt')


plt.plot(training_performance[:, 0], training_performance[:, 1])
plt.plot(training_performance_nokb[:, 0], training_performance_nokb[:, 1])
plt.title('Training Performance digits 3, 8, 5, 6')
plt.legend(['accuracy with kb', 'accuracy without kb'])
plt.show()


plt.plot(training_performance[:, 0], training_performance[:, 2])
plt.plot(training_performance_nokb[:, 0], training_performance_nokb[:, 2])
plt.title('Training Error digits 3, 8, 5, 6')
plt.legend(['error with kb', 'error without kb'])
plt.show()

plt.plot(training_performance_nokb_2[:, 0], training_performance_nokb_2[:, 1])
plt.title('Training Performance digits 3, 8')
plt.legend(['accuracy without kb'])
plt.show()

plt.plot(training_performance_nokb_2[:, 0], training_performance_nokb_2[:, 2])
plt.title('Training Performance digits 3, 8')
plt.legend(['error without kb'])
plt.show()


print('done')