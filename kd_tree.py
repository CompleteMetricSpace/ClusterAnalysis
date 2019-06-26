import random, math, numpy
from collections import deque


points_node_dict = {}
def get_all_points(node):
    if node in points_node_dict:
        return points_node_dict[node]
    result = []
    queue = deque([node])
    while queue:
        n = queue.pop()
        if n.is_leaf:
            result.extend(n.point_list)
        else:
            queue.extend([n.left_child, n.right_child])
    points_node_dict[node] = result
    return result

class KDTree:
    """
    KD-Tree class
    See wikipedia for more information
    """
    def __init__(self, root, sq_metric):
        self.root = root
        self.sq_metric = sq_metric


    def get_points_in_range(self, x, range):
        result = []
        self.root.get_points_in_range(x, range, result)
        return result


    def get_points_in_range_it(self, x, range_radius):
        """
        Determines the list of all points that are in the closed ball around x of radius range_radius
        :param x: a point (query point)
        :param range_radius: a positive number (query radius)
        :return: a list of all points in B(x, range_radius)
        """
        result = []
        queue = deque([self.root])
        while queue:
            node = queue.pop()
            if node.is_leaf:
                result.extend(node.get_points_in_range(x, range_radius, self.sq_metric))
            else:
                # First check if this node is completely included in the range ball around x
                # This is very inefficient in high dimensions

                #if is_contained_2_exact(x, range_radius, node.min_values, node.max_values):
                #    result.extend(get_all_points(node))
                #    continue
                x_val = x[node.axis]
                if abs(x_val - node.split_value) <= range_radius:
                    queue.append(node.left_child)
                    queue.append(node.right_child)
                else:
                    if x_val < node.split_value:
                        queue.append(node.left_child)
                    else:
                        queue.append(node.right_child)

        return result


    def construct_m_tree(self, m_tree):
        self.root.construct_m_tree(m_tree)


class InnerNode:

    def __init__(self, dim, left_child, right_child, split_value, axis, min_values, max_values):
        self.dim = dim
        self.left_child = left_child
        self.right_child = right_child
        self.split_value = split_value
        self.axis = axis
        self.min_values = min_values
        self.max_values = max_values
        self.is_leaf = False


    def get_points_in_range(self, x, range, sq_metric):
        x_val = x[self.axis]
        if abs(x_val-self.split_value) < range:
            result = self.left_child.get_points_in_range(x, range, sq_metric)
            result.extend(self.right_child.get_points_in_range(x, range, sq_metric))
            return result
        else:
            if x_val < self.split_value:
                return self.left_child.get_points_in_range(x, range, sq_metric)
            else:
                return self.right_child.get_points_in_range(x, range, sq_metric)


    def append_all_points(self, result):
        self.left_child.append_all_points(result)
        self.right_child.append_all_points(result)

    def construct_m_tree(self, m_tree):
        self.right_child.construct_m_tree(m_tree)
        self.left_child.construct_m_tree(m_tree)

class LeafNode:

    def __init__(self, dim, point_list, min_values, max_values):
        self.dim = dim
        self.point_list = point_list
        self.min_values = min_values
        self.max_values = max_values
        self.is_leaf = True


    def get_points_in_range(self, x, range_radius, sq_metric):
        r_sq = range_radius*range_radius
        return [y for y in self.point_list if sq_metric(x, y) <= r_sq]


    def append_all_points(self, result):
        result.extend(self.point_list)

    def construct_m_tree(self, m_tree):
        for p in self.point_list:
            m_tree.insert(p)


def construct_kd_tree_node(data_list, k, dim, max_number_of_points_in_leafs, min_values, max_values, sliding=False):
    if data_list.__len__() < max_number_of_points_in_leafs:
        return LeafNode(dim, data_list, min_values, max_values)
    else:
        if sliding:
            midpoint = (min_values[k]+max_values[k])/2
            left_data = []
            right_data = []
            for x in data_list:
                if x[k] < midpoint:
                    left_data.append(x)
                else:
                    right_data.append(x)
            if right_data == []:
                midpoint = max([x[k] for x in left_data])
            if left_data == []:
                midpoint = min([x[k] for x in right_data])

            left_min_values = min_values.copy()
            left_max_values = max_values.copy()
            left_max_values[k] = midpoint
            right_min_values = min_values.copy()
            right_max_values = max_values.copy()
            right_min_values[k] = midpoint

            left_node = construct_kd_tree_node(left_data, (k + 1) % dim, dim, max_number_of_points_in_leafs,
                                               left_min_values, left_max_values)
            right_node = construct_kd_tree_node(right_data, (k + 1) % dim, dim, max_number_of_points_in_leafs,
                                                right_min_values, right_max_values)
            return InnerNode(dim, left_node, right_node, midpoint, k, min_values, max_values)

        else:
            array = random.sample(data_list, min(30, max_number_of_points_in_leafs))
            median = numpy.median(list(map(lambda x: x[k], array)))
            left_data = []
            right_data = []
            for x in data_list:
                if x[k] < median:
                    left_data.append(x)
                else:
                    right_data.append(x)


            left_min_values = min_values.copy()
            left_max_values = max_values.copy()
            left_max_values[k] = median
            right_min_values = min_values.copy()
            right_max_values = max_values.copy()
            right_min_values[k] = median

            left_node = construct_kd_tree_node(left_data, (k+1) % dim, dim, max_number_of_points_in_leafs, left_min_values, left_max_values)
            right_node =  construct_kd_tree_node(right_data, (k+1) % dim, dim, max_number_of_points_in_leafs, right_min_values, right_max_values)
            return InnerNode(dim, left_node, right_node, median, k, min_values, max_values)




def create_kd_tree(point_list, leaf_points, sq_metric):
    dim = point_list[0].__len__()
    return KDTree(construct_kd_tree_node(point_list, 0, dim, leaf_points, [0 for i in range(dim)], [1 for i in range(dim)], sliding=True), sq_metric)
