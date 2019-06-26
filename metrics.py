import math

#Ok, this is quite embarassing, but we have to do this
#since python function calls are really slow

dataList = None

dist_dict = {}
def dist_mem(x, y):
    global reused
    s = frozenset([x,y])
    if s in dist_dict:
        return dist_dict[s]
    else:
        d = math.sqrt(sum([(a-b)*(a-b) for a,b in zip(x,y)]))
        dist_dict[s] = d
        return d


def sq_dist(x, y):
    return sum([(a-b)*(a-b) for a,b in zip(x,y)])

def sq_dist_2(x, y):
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])

def sq_dist_3(x, y):
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])

def sq_dist_4(x, y):
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])+(x[3]-y[3])*(x[3]-y[3])

def sq_dist_9(x, y):
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])\
           +(x[3]-y[3])*(x[3]-y[3])+(x[4]-y[4])*(x[4]-y[4])+(x[5]-y[5])*(x[5]-y[5]) \
           + (x[6] - y[6]) * (x[6] - y[6]) + (x[7] - y[7]) * (x[7] - y[7]) + (x[8] - y[8]) * (x[8] - y[8])

def sq_dist_10(x, y):
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])\
           +(x[3]-y[3])*(x[3]-y[3])+(x[4]-y[4])*(x[4]-y[4])+(x[5]-y[5])*(x[5]-y[5]) \
           + (x[6] - y[6]) * (x[6] - y[6]) + (x[7] - y[7]) * (x[7] - y[7]) + (x[8] - y[8]) * (x[8] - y[8]) \
           + (x[9]-y[9])*(x[9]-y[9])


def dist_mem_index(i,j):
    return dist_mem(dataList[i], dataList[j])

def sq_dist_index(i, j):
    return sum([(a-b)*(a-b) for a,b in zip(dataList[i],dataList[j])])

def sq_dist_2_index(i, j):
    x, y = dataList[i], dataList[j]
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])

def sq_dist_3_index(i, j):
    x, y = dataList[i], dataList[j]
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])

def sq_dist_4_index(i, j):
    x, y = dataList[i], dataList[j]
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])+(x[3]-y[3])*(x[3]-y[3])

def sq_dist_10_index(i, j):
    x, y = dataList[i], dataList[j]
    return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])\
           +(x[3]-y[3])*(x[3]-y[3])+(x[4]-y[4])*(x[4]-y[4])+(x[5]-y[5])*(x[5]-y[5]) \
           + (x[6] - y[6]) * (x[6] - y[6]) + (x[7] - y[7]) * (x[7] - y[7]) + (x[8] - y[8]) * (x[8] - y[8]) \
           + (x[9]-y[9])*(x[9]-y[9])



def choose_sq_metric(dimension):
    if dimension == 2:
        return sq_dist_2
    if dimension == 3:
        return sq_dist_3
    if dimension == 4:
        return sq_dist_4
    if dimension == 9:
        return sq_dist_9
    if dimension == 10:
        return sq_dist_10
    return sq_dist


def choose_sq_metric_index(dimension):
    if dimension == 2:
        return sq_dist_2_index
    if dimension == 3:
        return sq_dist_3_index
    if dimension == 4:
        return sq_dist_4_index
    if dimension == 10:
        return sq_dist_10_index
    return sq_dist_index


def choose_metric(dimension):
    return dist_mem


def choose_metric_index(dimension):
    return dist_mem_index