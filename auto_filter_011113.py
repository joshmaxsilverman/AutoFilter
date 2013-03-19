import os
import string
import sets
import csv
import re
# import numpypy
import numpy
import math
import sys
import time
# import matplotlib.pyplot as plt

path = '/Users/joshsilverman/Dropbox/Research/AutoFilter/gogat40c'

os.chdir(path)

# 1 - Get list of file stems
fileList = os.listdir(path)[2:]
stems = list(set(map( lambda file: re.split(".[a-z]", file)[0], fileList)))

# 2 - Function to get relevant fit and data file information
def setup(stem):
    # Get the fit point pairs
    fit_file = open(stem + ".fit")
    fit = numpy.array([map(float, item) for item in csv.reader(fit_file)])
    # Get the data point pairs    
    dat_file = open(stem + ".dat")
    data = numpy.array([map(float, item) for item in csv.reader(dat_file)])
    
    return fit, data

# 3 - Function to mask the data outside the fit and return a residual score
def masking(stem):

    time2 = time.time()
    
    fit, data = setup(stem)
    length_data, length_fit = len(data), len(fit)    
    
    # 3a - Get the length of the fit and take 2%
    fit_tail_length = int(len(fit) * 0.02)
    # 3b - Get the intensities for the tail piece
    fit_tail = fit[-fit_tail_length:-1].transpose()[1]
    # 3c - Calculate the median value for the tail region
    fit_baseline = numpy.mean(fit_tail)
    # 3d - Calculate the standard deviation of the tail region
    fit_noise = math.sqrt(numpy.var(fit_tail))
    
    checkpoints = []
    # 3e - If fit point is far above baseline, return data index
    for index in range(length_data):
        if cmp(fit[index * int(length_fit / length_data)][1] - (fit_baseline + 150 * fit_noise), 0) == 1:
            checkpoints.append(index)
        else:
            pass

    # 3f - Scale the checked indices into fit indices
    checked_scaled = [item * int(length_fit / length_data) for item in checkpoints]
    # 3g - Return the relevant unmasked fit and data points
    relevant_data_points = numpy.array([data[item][1] for item in checkpoints])
    relevant_fit_points = numpy.array([fit[item][1] for item in checked_scaled])
    # 3h - Calculate the difference, square the sum and divide the magnitude by total area under fit curve
    deltas = relevant_data_points - relevant_fit_points
    residue = sum( delta**2 for delta in deltas)
    score = math.sqrt(residue) / sum(relevant_fit_points)    
    time1 = time.time()

    print "Process took %0.2f seconds, score was %0.2f" % (time1 - time2, score)
    
#     return stem, score
    return time.time() - time1

times = map(lambda item: masking(item), stems[0:50])
print numpy.mean(times)

# To check the unmasked points

# plt.plot(checked_scaled, 10 * numpy.ones(len(checked_scaled)), 'ko', ms = 2)
# plt.plot(range(length_fit), fit.transpose()[1])
# plt.show()
        
