import numpy as np
from PIL import Image, ImageDraw

img_width = 28
visualization_file_path: str = "../cmake-build-debug/remote/archived/output-3/4digits_1/visualization.txt"
cl_file_path = "../cmake-build-debug/remote/archived/output-3/4digits_1/classifier.txt"
cf_file_path = "../cmake-build-debug/remote/archived/output-3/4digits_1/code_fragment.txt"
filter_file_path = "../cmake-build-debug/remote/archived/output-3/4digits_1/filter.txt"

# load classifiers ids and their code fragment ids
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


def get_image(val):
    data = np.zeros((28, 28))
    data += val
    return data


img = get_image(0)
img_l = get_image(1)
img_u = get_image(0)


def load_filter(id, size):
    f_id = -1
    f = open(filter_file_path)
    line = f.readline()
    while line:
        tokens = line.strip().split()
        f_id = int(tokens[1])
        if f_id == id:
            break
        line = f.readline()
        line = f.readline()
        line = f.readline()

    assert f_id != -1
    line = f.readline()
    tokens = line.strip().split()
    lb = []
    ub = []
    for i in range(size*size+1):
        if i == 0:
            continue
        lb.append(float(tokens[i]))
    line = f.readline()
    tokens = line.strip().split()
    for i in range(size*size+1):
        if i == 0:
            continue
        ub.append(float(tokens[i]))

    return lb, ub


# update lower and upper bounds from filter
def update_bounds(id, start_x, start_y, size, dilated):
    lb, ub = load_filter(id, size)
    step = 1
    if dilated:
        step = 2
        effective_size = size + size - 1

    for y in range(size):
        for x in range(size):
            if img_l[(start_y+y*step), start_x+x*step] > lb[y*size + x]:
                img_l[(start_y+y*step), start_x+x*step] = lb[y*size + x]
            if img_u[(start_y+y*step), start_x+x*step] < ub[y*size + x]:
                img_u[(start_y+y*step), start_x+x*step] = ub[y*size + x]


def get_pixel_color(x,y):
    color = (255, 0, 0)
    if img_l[x, y] == 1:  # if the pixel interval has not be initialized then its don't care
        return color
    if img_l[x, y] < 0.25 and img_u[x, y] > 0.75:  # wide interval means don't care
        return color

    mid = (img_l[x,y] + img_u[x,y]) / 2
    # real to 255 scale
    c = int(mid*255)
    color = (c, c, c)

    # color = "ff0000"
    # if mid < 0.25:
    #     color = "000000"  # black
    # elif mid < 0.5:
    #     color = "D3D3D3"  # light grey

    return color


def visualize_intervals(dc):
    for y in range(img_l.shape[0]):
        for x in range(img_l.shape[0]):
            dc.point((x, y), get_pixel_color(x,y))


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
                        # top right corner
                        y = position // img_width
                        x = position % img_width
                        update_bounds(filter, x, y, size, dilated)

    # visualize now
    base_img = Image.fromarray(img).convert("RGB")
    dc = ImageDraw.Draw(base_img)  # draw context
    visualize_intervals(dc)

    print("filters drawn: "+str(filters_drawn))
    base_img = base_img.resize((300, 300))
    base_img.show()


for i in range(100):
    visualize_image(i, True)
    input("press any key to continue")





print('done')
