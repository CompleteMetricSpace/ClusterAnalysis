import time, csv, random, math
from collections import deque
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import metrics as mtr



class ClusterList:

    def __init__(self, root, metric_sq):
        self.root = root
        self.metric_sq = metric_sq


    def search(self, obj, range_radius):
        range_radius_sq = range_radius*range_radius
        result = []
        current_node = self.root
        while current_node != None:
            dist_sq = self.metric_sq(obj, current_node.center)
            if dist_sq > (range_radius+current_node.radius)*(range_radius+current_node.radius):
                current_node = current_node.child_node
                continue
            if dist_sq <= range_radius_sq:
                result.append(current_node.center)
            if range_radius >= current_node.radius and dist_sq <= (range_radius-current_node.radius)*(range_radius-current_node.radius):
                result.extend([x[0] for x in current_node.element_list])
            else:
                for (y, d_sq) in current_node.element_list:
                    if range_radius_sq >= dist_sq + d_sq and 4*dist_sq*d_sq <= (range_radius_sq - dist_sq - d_sq)*(range_radius_sq - dist_sq - d_sq):
                        result.append(y)
                        continue
                    r_sq = self.metric_sq(obj, y)
                    if r_sq <=range_radius_sq:
                        result.append(y)
            if current_node.radius >= range_radius and dist_sq <= (current_node.radius - range_radius)*(current_node.radius - range_radius):
                break #Range disk was completely included in this node
            current_node = current_node.child_node
        return result


    def get_average_elem_number(self):
        sum = 0
        counter = 0
        current_node = self.root
        while not current_node is None:
            sum += current_node.element_list.__len__()
            counter += 1
            current_node = current_node.child_node
        return (sum/counter, sum, counter)

class ClusterNode:

    def __init__(self, child_node, radius, center, element_list):
        self.child_node = child_node
        self.radius = radius
        self.center = center
        self.element_list = element_list




def build_cluster_node(point_list, radius, metric_sq):
    if point_list == []:
        return None
    c = random.sample(point_list,1)[0]
    element_list = []
    new_point_list = []
    for p in point_list:
        if p == c:
            continue
        d_sq = metric_sq(p,c)
        if d_sq <= radius*radius:
            element_list.append((p, d_sq))
        else:
            new_point_list.append(p)
    return ClusterNode(build_cluster_node(new_point_list, radius, metric_sq), radius, c, element_list)


def build_cluster_node_it(point_list, radius, metric_sq):
    new_point_list = {p:True for p in point_list}
    c = random.sample(list(new_point_list), 1)[0]
    del new_point_list[c]
    element_list = []
    tmp = {}
    for p in new_point_list:
        d_sq = metric_sq(p, c)
        if d_sq <= radius * radius:
            element_list.append((p, d_sq))
        else:
            tmp[p] = True
    new_point_list = tmp
    current_node = ClusterNode(None, radius, c, element_list)
    root_node = current_node
    while new_point_list.__len__() > 0:
        c = random.sample(list(new_point_list), 1)[0]
        del new_point_list[c]
        element_list = []
        tmp = {}
        for p in new_point_list:
            d_sq = metric_sq(p, c)
            if d_sq <= radius * radius:
                element_list.append((p, d_sq))
            else:
                tmp[p] = True
        new_point_list = tmp
        current_node.child_node = ClusterNode(None, radius, c, element_list)
        current_node = current_node.child_node
    return root_node


def build_cluster_list(point_list, radius, metric_sq):
    return ClusterList(build_cluster_node(point_list, radius, metric_sq), metric_sq)


def build_cluster_list_it(point_list, radius, metric_sq):
    return ClusterList(build_cluster_node_it(point_list, radius, metric_sq), metric_sq)