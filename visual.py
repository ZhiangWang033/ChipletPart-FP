import os
import matplotlib.pyplot as plt
from math import log




#sample 3
#net_file = "./test/4_14_4_400_1/chiplet_netlist_64.hgr"

net_file = "./test/example_3/sample.hgr"


shape_file = net_file + ".output"




with open(shape_file) as f:
    lines = f.readlines()
f.close()

lx_list = []
ly_list = []
width_list = []
height_list = []

for line in lines:
    words = line.split()
    lx_list.append(float(words[0]))
    ly_list.append(float(words[1]))
    width_list.append(float(words[2]))
    height_list.append(float(words[3]))

net_src_list = []
net_dst_list = []
net_weight_list = []
net_reach_list = []

with open(net_file) as f:
    lines = f.readlines()
f.close()

num_nets = int(lines[0].split()[0])
for i in range(1, num_nets + 1):
    words = lines[i].split()
    net_src_list.append(int(words[0]))
    net_dst_list.append(int(words[1]))
    net_weight_list.append(int(words[2]))
    net_reach_list.append(int(words[3]))

plt.figure()
for i in range(len(lx_list)):
    lx = lx_list[i]
    ly = ly_list[i]
    width = width_list[i]
    height = height_list[i]
    rectangle = plt.Rectangle((lx, ly), width, height, fc = "green", ec = "blue")
    plt.gca().add_patch(rectangle)

plt.axis('equal')
plt.show()














