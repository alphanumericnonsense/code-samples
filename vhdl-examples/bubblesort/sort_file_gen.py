import random
import math
import sys

depth_g = 10
width_g = 8
num_trials_g = 10

if len(sys.argv) == 4:
    try:
        depth_g = int(sys.argv[1])
        width_g = int(sys.argv[2])
        num_trials_g = int(sys.argv[3])
    except:
        print("something wrong with cmdline args, using defaults.")


print(f"creating sort.txt with {num_trials_g}*{depth_g} {width_g}-bit integers, one on each line...")

with open("sort.txt", 'w') as f:
    for i in range(depth_g*num_trials_g):
        x = math.floor((2**width_g)*random.random())
        f.write(str(x))
        f.write("\n")

print("... finished.")
