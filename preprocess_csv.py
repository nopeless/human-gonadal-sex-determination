import csv
import numpy as np

# np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0], np.int32).tofile("groups.bin")
# exit(0)

reader = csv.reader(open("matlab.csv", "r"))

data = []

next(reader)

for idx, line in enumerate(reader):
    numbers = [int(x) - 1 for x in line]

    if max(numbers) != 1:
        continue

    # print(numbers)

    data.append(numbers)

# data = data[:1]

np.array(data, np.int32).tofile("groups.bin")
