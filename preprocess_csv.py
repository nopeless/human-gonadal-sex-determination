import csv
import numpy as np

reader = csv.reader(open("matlab.csv", "r"))

data = []

next(reader)

for idx, line in enumerate(reader):
    numbers = [int(x) - 1 for x in line]

    # print(numbers)

    data.append(numbers)

np.array(data, np.int32).tofile("groups.bin")
