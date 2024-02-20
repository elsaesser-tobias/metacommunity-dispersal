import nlmpy
import np
import random
import csv
import skimage
import math

stream = 1
ag = 0.005
forest = 0.00017
urban = 0.0000032

#import seed number from R script and set seed for nlmpy
with open('py_seed.txt') as f:
    seed = f.readlines()
seed = int(seed[0])
np.random.seed(seed)

with open('py_cellsz.txt') as f:
    sz = f.readlines()
sz = int(sz[0])


#import habitat x any y coordinates from R script
with open('x_points.txt') as f:
    x_points = f.readlines()
x_points = list(map(int, x_points))

with open('y_points.txt') as f:
    y_points = f.readlines()
y_points = list(map(int, y_points))

#import stream order and downnode
with open('SO.txt') as f:
    SO = f.readlines()
SO = list(map(int, SO))


with open('DN.txt') as f:
    DN = f.readlines()
DN = list(map(int, DN))

#import habitat points
with open('habitatsx.txt') as f:
    habitatsx = f.readlines()
habitatsx = list(map(int, habitatsx))

with open('habitatsy.txt') as f:
    habitatsy = f.readlines()
habitatsy = list(map(int, habitatsy))

with open('urban_landscape_percent.txt') as f:
    urban_percent = f.readlines()
urban_percent = float(urban_percent[0])

with open('forest_landscape_percent.txt') as f:
    forest_percent = f.readlines()
forest_percent = float(forest_percent[0])

with open('agriculture_landscape_percent.txt') as f:
    ag_percent = f.readlines()
ag_percent = float(ag_percent[0])

# create a mask of different habitat qualitites, with the randomElementNN function
#nlmElement = nlmpy.randomElementNN(10000, 10000, 500)
nlmElement_X = nlmpy.randomClusterNN(500, 500, 0.3825, n='8-neighbourhood')
nlmElement = np.repeat(np.repeat(nlmElement_X, 20, axis = 0), 20, axis = 1)
#nlmElement = nlmpy.random(10000, 10000)

for i in range(0,nlmElement.shape[0]):
    for j in range(0,nlmElement.shape[1]):
        ele = nlmElement[i,j]
        if ele <= ag_percent:
            nlmElement[i,j] = ag
        if ele <= forest_percent + ag_percent and ele > ag_percent:
            nlmElement[i,j] = forest
        if ele <= forest_percent + ag_percent + urban_percent and ele > forest_percent + ag_percent:
            nlmElement[i,j] = urban


#connect streams, use stream order for how broad
nlmElement[x_points[0] - 1,y_points[0] - 1] = stream
check = 0
for i in range(0,len(DN)):
    if(DN[i] == 0):
        continue
    #print("I:",i)
    nlmElement[x_points[i] - 1,y_points[i] -1] = stream
    # check orientation of downstreamnode
    x_step = 1
    y_step = 1
    if(x_points[i] - x_points[DN[i] -1] < 0):
        x_step = -1
    if(y_points[i] - y_points[DN[i] -1] < 0):
        y_step = -1   
    x_range = range(0,x_points[i] - x_points[DN[i] -1], x_step)
    y_range = range(0,y_points[i] - y_points[DN[i] -1], y_step)
    #if either x or y range is zero, create zero vector for the distance in that direction
    if(len(x_range) < 100):
       x_range = np.zeros(100, dtype = int)
    if(len(y_range) < 100):
       y_range = np.zeros(100, dtype = int)

    for j in range(0,len(x_range)):
        nlmElement[x_points[i] - x_range[j] - 1,y_points[i] - y_range[j] - 1] = stream
        for k in range(- SO[i] * 5,SO[i] * 5):
            nlmElement[x_points[i] - x_range[j] + k - 1,y_points[i] - y_range[j] - 1] = stream
            for m in range(- SO[i] * 5,SO[i] * 5):
                nlmElement[x_points[i] - x_range[j] + k - 1,y_points[i] - y_range[j] + m - 1] = stream
            


# Specify the CSV file path
csv_file_path = 'nlmElement.csv'

# Write the matrix to a CSV file
with open(csv_file_path, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerows(nlmElement)


