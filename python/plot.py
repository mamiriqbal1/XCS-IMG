import numpy as np
import matplotlib.pyplot as plt


training_performance_nokb_2 = np.loadtxt('../cmake-build-debug/output_2_6_2000_nokb_save/output_training.txt')


plt.plot(training_performance_nokb_2[:, 0], training_performance_nokb_2[:, 1])
plt.title('Training Accuracy digits 3, 8')
plt.legend(['Training Accuracy digits 3, 8'])
plt.show()

plt.plot(training_performance_nokb_2[:, 0], training_performance_nokb_2[:, 2])
plt.title('Training Error digits 3, 8')
plt.legend(['Training Error digits 3, 8'])
plt.show()


print('done')