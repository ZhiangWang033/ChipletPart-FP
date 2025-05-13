import os
import matplotlib.pyplot as plt


file_name = "reach_sweep_chiplet_25.txt"

with open(file_name, "r") as f:
    content = f.read().splitlines()
f.close()

reach = []
package_size = []

for line in content:
    reach.append(int(line.split(" ")[0]))
    package_size.append(int(line.split(" ")[1]))
    
plt.plot(reach, package_size, "*")
plt.xlabel("Reach")
plt.ylabel("Package Size")
plt.title("Reach vs Package Size (25 Chiplets)")
plt.show()


