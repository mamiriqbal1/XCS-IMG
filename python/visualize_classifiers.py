import numpy as np
from PIL import Image, ImageDraw

img_width = 28
visualization_file_path = "../cmake-build-debug/remote/output/4digits_1/visualization.txt"
image_file_path = "../data/mnist/mnist_test_3_8_5_6.txt"
cl_file_path = "../cmake-build-debug/remote/output/4digits_1/classifier.txt"
cf_file_path = "../cmake-build-debug/remote/output/4digits_1/code_fragment.txt"

# load classifiers is and their code fragment ids
cl_cf = {}
f = open(cl_file_path)
line = f.readline()
while line:
    tokens = line.strip().split()
    cl_id = int(tokens[1])
    line = f.readline()
    tokens = line.strip().split()
    cf = []
    for token in tokens:
        cf.append(int(token))
    cl_cf[cl_id] = cf
    line = f.readline()
f.close()

# load code fragments and their filter ids
cf_filter = {}
f = open(cf_file_path)
line = f.readline()
while line:
    tokens = line.strip().split()
    cl_id = int(tokens[0])
    filters = []
    for token in tokens:
        if token.startswith("D"):
            filters.append(int(token[1:]))
    cf_filter[cl_id] = filters
    line = f.readline()
f.close()


def get_image(img_id):
    img_file = np.loadtxt(image_file_path)
    item = img_file[img_id]
    img_class = int(item[-1])
    data = item[:-1]
    # denormalize
    data = data * 255
    data = data.reshape(28, 28)
    return img_class, data


cl_clclass = []
filter_position = {}


def visualize_image(img_id, rectangle):
    # load visualization data
    actual_class = -1
    predicted_class = -1
    f = open(visualization_file_path)
    line = f.readline()
    while line:
        # skip img_id lines to get to the the right image
        for i in range(img_id):
            line = f.readline()
            line = f.readline()
            line = f.readline()
        tokens = line.strip().split()
        read_img_id = int(tokens[0])
        assert(read_img_id == img_id)
        actual_class = int(tokens[1])
        predicted_class = int(tokens[2])
        print("actual class: " + str(actual_class) + " predicted class: " + str(predicted_class) + "\n")
        line = f.readline()  # classifier_id predicted_class ...
        tokens = line.strip().split()
        i = 0
        while i < len(tokens):
            classifier = tokens[i]
            i += 1
            classifier_class = tokens[i]
            i += 1
            cl_clclass.append((int(classifier), int(classifier_class)))
        line = f.readline()  # filter_id matched_position
        tokens = line.strip().split()
        i = 0
        while i < len(tokens):
            filter_id = int(tokens[i])
            i += 1
            position = int(tokens[i])
            i += 1
            size = int(tokens[i])
            i += 1
            dilated = int(tokens[i])
            i += 1
            filter_position[filter_id] = (position, size, dilated)
        break
    f.close()
    img = get_image(img_id)
    assert(img[0] == actual_class)
    base_img = Image.fromarray(img[1]).convert("RGB")
    dc = ImageDraw.Draw(base_img)  # draw context
    # draw dots/rectangles based on filter position from positive classifiers
    filters_drawn = 0
    for classifier in cl_clclass:
        if classifier[1] == actual_class:  # if positive classifier
            classifier_id = classifier[0]
            # get classifier code fragments
            code_fragments = cl_cf[classifier_id]
            for cf in code_fragments:
                filters = cf_filter[cf]  # filter
                for filter in filters:
                    position = filter_position[filter][0]
                    size = filter_position[filter][1]
                    dilated = filter_position[filter][2]
                    if position >= 0:
                        filters_drawn += 1
                        if dilated:
                            size += (size-1)
                        # top right corner
                        y = position // img_width
                        x = position % img_width
                        if rectangle:
                            shape = [(x,y), (x+size,y+size)]
                            dc.rectangle(shape, outline="#ff0000")
                        else:
                            # center point
                            x += size//2
                            y += size//2
                            dc.point((x, y), fill="#ff0000")
    print("filters drawn: "+str(filters_drawn))
    base_img = base_img.resize((300, 300))
    base_img.show()


for i in range(100):
    visualize_image(i, True)
    input("press any key to continue")





print('done')
