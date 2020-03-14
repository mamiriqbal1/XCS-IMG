from keras.datasets import mnist
import numpy as np


def prepare_mnist_all(out_path):
    (x_train, y_train), (x_test, y_test) = mnist.load_data()
    # reshape
    x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
    x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))
    # normalize between 0 and 1
    # x_train = (x_train - x_train.min()) / (x_train.max() - x_train.min())
    # x_test = (x_test - x_test.min()) / (x_test.max() - x_test.min())

#    np.savetxt(out_path + 'all_train_x_mnist.txt', x_train)
#    np.savetxt(out_path + 'all_test_x_mnist.txt', x_test)
#
#    np.savetxt(out_path + 'all_train_y_mnist.txt', y_train, fmt='%d')
#    np.savetxt(out_path + 'all_test_y_mnist.txt', y_test, fmt='%d')

    np.savetxt(out_path + 'all_train_mnist.txt', np.hstack((x_train, np.reshape(y_train, (-1, 1)))))
    np.savetxt(out_path + 'all_test_mnist.txt', np.hstack((x_test, np.reshape(y_test, (-1, 1)))))

def prepare_mnist_3_8(out_path):
    (x_train, y_train), (x_test, y_test) = mnist.load_data()

    # prepare the train set
    # extract sample for 3 & 8
    index_3 = np.where(y_train == 3)[0]
    index_8 = np.where(y_train == 8)[0]
    index_trim = np.concatenate((index_3, index_8), axis=None)
    index_trim = np.sort(index_trim)
    y_train = y_train[index_trim]
    x_train = x_train[index_trim]

    x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))

    # we are working with 3 and 8 only, we map 3 to 0 and 8 to 1
    y_train[y_train == 3] = 0
    y_train[y_train == 8] = 1

    # prepare the test set
    # extract sample for 3 & 8
    index_3 = np.where(y_test == 3)[0]
    index_8 = np.where(y_test == 8)[0]
    index_trim = np.concatenate((index_3, index_8), axis=None)
    index_trim = np.sort(index_trim)
    y_test = y_test[index_trim]
    x_test = x_test[index_trim]

    x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))

    # we are working with 3 and 8 only, we map 3 to 0 and 8 to 1
    y_test[y_test == 3] = 0
    y_test[y_test == 8] = 1

    # # normalize data
    maxv = max(x_train.max(), x_test.max())
    minv = min(x_train.min(), x_test.min())
    x_train = (x_train - minv) / (maxv - minv)
    x_test = (x_test - minv) / (maxv - minv)

    np.savetxt(out_path + 'mnist_train_3_8.txt', np.hstack((x_train, np.reshape(y_train, (-1, 1)))))
    np.savetxt(out_path + 'mnist_test_3_8.txt', np.hstack((x_test, np.reshape(y_test, (-1, 1)))))


def prepare_mnist_digits(out_path, digits, postfix):
    (x_train, y_train), (x_test, y_test) = mnist.load_data()

    # prepare the train set
    index_all = np.ndarray(shape=(0,), dtype=np.int64)
    for digit in digits:
        index_digit = np.where(y_train == digit)[0]
        index_all = np.concatenate((index_all, index_digit), axis=None)

    index_all = np.sort(index_all)
    y_train = y_train[index_all]
    x_train = x_train[index_all]

    x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))

    for digit, label in zip(digits, range(len(digits))):
        y_train[y_train == digit] = label

    # prepare the test set
    index_all = np.ndarray(shape=(0,), dtype=np.int64)
    for digit in digits:
        index_digit = np.where(y_test == digit)[0]
        index_all = np.concatenate((index_all, index_digit), axis=None)

    index_all = np.sort(index_all)
    y_test = y_test[index_all]
    x_test = x_test[index_all]

    x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))

    for digit, label in zip(digits, range(len(digits))):
        y_test[y_test == digit] = label


    # # normalize data
    maxv = max(x_train.max(), x_test.max())
    minv = min(x_train.min(), x_test.min())
    x_train = (x_train - minv) / (maxv - minv)
    x_test = (x_test - minv) / (maxv - minv)

    np.savetxt(out_path + 'train_mnist' + postfix + '.txt', np.hstack((x_train, np.reshape(y_train, (-1, 1)))))
    np.savetxt(out_path + 'test_mnist' + postfix + '.txt', np.hstack((x_test, np.reshape(y_test, (-1, 1)))))


prepare_mnist_digits('../data/mnist/', [3, 8], '_3_8')
prepare_mnist_digits('../data/mnist/', [3, 8, 5, 6], '_3_8_5_6')

#prepare_mnist_3_8('../data/mnist/')

#repare_mnist_all('../data/mnist/')
