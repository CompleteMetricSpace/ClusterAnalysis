from operator import mul
from functools import reduce
from collections import deque
import random
from morton import morton_code
import math


node_point_dict = {}
def get_points(node):
    if node in node_point_dict:
        return node_point_dict[node]
    points = []
    queue = deque([node])
    while queue:
        n = queue.pop()
        if n.is_leaf:
            points.extend(n.entry_list)
        else:
            queue.extend(n.subtree)
    node_point_dict[node] = points
    return points



class Rect:

    def __init__(self, min_dimensions, max_dimensions):
        self.min_dimensions = min_dimensions
        self.max_dimensions = max_dimensions
        self.dim = min_dimensions.__len__()


    def get_dim(self):
        return self.dim

    def get_center(self):
        return tuple((a+b)/2 for a,b in zip(self.min_dimensions, self.max_dimensions))

    def get_volume(self):
        return reduce(mul, [b - a for a, b in list(zip(self.min_dimensions, self.max_dimensions))], 1)

    def contains(self, obj):
        return all((obj[i] >= self.min_dimensions[i] and obj[i] <= self.max_dimensions[i] for i in range(obj.__len__())))

    def contains_rect(self, rect):
        return all((rect.min_dimensions[i] >= self.min_dimensions[i] and rect.max_dimensions[i] <= self.max_dimensions[i]  for i in range(rect.get_dim())))

    def does_intersect(self, rect):
        return not any((self.max_dimensions[i] < rect.min_dimensions[i] or self.min_dimensions[i] > rect.max_dimensions[i] for i in range(self.get_dim())))


def get_intersection(rect_1, rect_2):
    return Rect(list(map(max, zip(rect_1.min_dimensions, rect_2.min_dimensions))),
                list(map(min, zip(rect_1.max_dimensions, rect_2.max_dimensions))))


def get_union(rect_1, rect_2):
    return Rect(list(map(min, zip(rect_1.min_dimensions, rect_2.min_dimensions))),
                list(map(max, zip(rect_1.max_dimensions, rect_2.max_dimensions))))


def get_union_list(rect_list):
    rect = rect_list[0]
    for r in rect_list[1:]:
        rect = get_union(rect, r)
    return rect


def get_union_points(obj_1, obj_2):
    return Rect(list(map(min, zip(obj_1, obj_2))), list(map(max, zip(obj_1, obj_2))))


def get_union_points_list(obj_list):
    rect = Rect(obj_list[0], obj_list[0])
    for o in obj_list[1:]:
        rect = get_enlarged_rect(rect, o)
    return rect

def get_enlarged_rect(rect, point):
    return Rect(list(map(min, zip(rect.min_dimensions, point))), list(map(max, zip(rect.max_dimensions, point))))


def get_overlap(rect_list):
    #Compute areas of intersection
    n = rect_list.__len__()
    area_dict = {}
    for i in range(n):
        for j in range(i):
            area_dict[frozenset([i,j])] = get_intersection(rect_list[i], rect_list[j]).get_volume()
    #Compute overlap for all nodes:
    overlap = []
    for i in range(n):
        overlap.append(sum([area_dict[frozenset([j,i])] for j in range(n) if j != i]))
    return overlap


def get_overlap_single(rect_list, index):
    n = rect_list.__len__()
    return sum([get_intersection(rect_list[index], rect_list[j]).get_volume() for j in range(n) if j != index])


def does_intersect(obj, radius, rect:Rect):
    #Optimized for dimension 2
    if rect.get_dim() == 2:
        return not (obj[0] - radius > rect.max_dimensions[0] or obj[0] + radius < rect.min_dimensions[0]
                    or obj[1] - radius > rect.max_dimensions[1] or obj[1] + radius < rect.min_dimensions[1])
    for i in range(obj.__len__()):
        if obj[i] - radius > rect.max_dimensions[i] or obj[i] + radius < rect.min_dimensions[i]:
            return False
    return True

def is_contained(obj, radius, rect):
    z = radius/math.sqrt(obj.__len__())
    return all([obj[i] + z >= rect.max_dimensions[i] and obj[i] - z <= rect.min_dimensions[i] for i in range(obj.__len__())])
  #  return (obj[0] + z >= rect.max_dimensions[0] and obj[0] - z <= rect.min_dimensions[0] and obj[1] + z >= rect.max_dimensions[1] and obj[1] - z <= rect.min_dimensions[1])


def get_all_points(node):
    result = []
    queue = deque([node])
    while queue:
        n = queue.pop()
        if n.is_leaf:
            result.extend(n.entry_list)
        else:
            queue.extend(n.subtree)
    return result

class RTree:
    """
    R-Tree class
    For more information see the original paper and the STR-paper.
    """
    def __init__(self, root, max_inner_node_children, min_inner_node_children, max_leaf_node_elements, min_leaf_node_elements, sq_metric):
        self.root = root
        self.max_inner_node_children = max_inner_node_children
        self.min_inner_node_children = min_inner_node_children
        self.max_leaf_node_elements = max_leaf_node_elements
        self.min_leaf_node_elements = min_leaf_node_elements
        self.sq_metric = sq_metric

    @staticmethod
    def get_next_node(node_list, rect_1, rect_2):
        vol_1 = rect_1.get_volume()
        vol_2 = rect_2.get_volume()
        volume_increase_difference = -math.inf
        next_node = None
        for node in node_list:
            d_1 = get_union(rect_1, node.rect).get_volume() - vol_1
            d_2 = get_union(rect_2, node.rect).get_volume() - vol_2
            d = abs(d_1-d_2)
            if volume_increase_difference < d:
                volume_increase_difference = d
                next_node = node
        return next_node


    @staticmethod
    def get_group_to_insert_node(node, rect_1, rect_2):
        vol_1 = rect_1.get_volume()
        vol_2 = rect_2.get_volume()
        inc_1 = get_union(node.rect, rect_1).get_volume()-vol_1
        inc_2 = get_union(node.rect, rect_2).get_volume()-vol_2
        if inc_1 < inc_2:
            return 1
        if inc_1 > inc_2:
            return 2
        #print("Tie on insert!")
        return 1

    @staticmethod
    def get_group_to_insert_obj(obj, rect_1, rect_2):
        vol_1 = rect_1.get_volume()
        vol_2 = rect_2.get_volume()
        inc_1 = get_enlarged_rect(rect_1, obj).get_volume() - vol_1
        inc_2 = get_enlarged_rect(rect_2, obj).get_volume() - vol_2
        if inc_1 < inc_2:
            return 1
        if inc_1 > inc_2:
            return 2
        #print("Tie on insert!")
        return 1

    @staticmethod
    def get_next_obj(obj_list, rect_1, rect_2):
        vol_1 = rect_1.get_volume()
        vol_2 = rect_2.get_volume()
        volume_increase_difference = -math.inf
        next_obj = None
        for obj in obj_list:
            d_1 = get_enlarged_rect(rect_1, obj).get_volume() - vol_1
            d_2 = get_enlarged_rect(rect_2, obj).get_volume() - vol_2
            d = abs(d_1-d_2)
            if volume_increase_difference < d:
                volume_increase_difference = d
                next_obj = obj
        return next_obj

    def insert(self, obj):
        if self.root == None:
            self.root = LeafNode(self, Rect(obj, obj), [obj])
        else:
            self.root.rect = get_enlarged_rect(self.root.rect, obj)
            split_node = self.root.insert(obj)
            if split_node != None:
                new_bounding_rect = get_union(self.root.rect, split_node.rect)
                self.root = InnerNode(self, new_bounding_rect, [self.root, split_node])


    def is_correct(self):
        return self.root.is_correct()


    def get_points_in_range_it(self, obj, range_radius):
        """
        This method computes all points within a radius of an object
        :param obj: a point
        :param range_radius: a positive number
        :return: the set of all points that are in the closed ball B(obj,range_radius)
        """
        result = []
        range_sq = range_radius*range_radius
        queue = deque([self.root])
        while queue:
            node = queue.pop()
            #Is contained is really inefficient in high dimensions
            #if is_contained(obj, range_radius, node.rect):
            #    result.extend(get_points(node))
            #    continue
            if node.is_leaf:
                result.extend([y for y in node.entry_list if self.sq_metric(y,obj) <= range_sq])
            else:
                queue.extend([n for n in node.subtree if does_intersect(obj, range_radius, n.rect)])
        return result

    def get_points_in_rect_it(self, rect:Rect):
        """
        This method computes all points within a radius of an object
        :param obj: a point
        :param range_radius: a positive number
        :return: the set of all points that are in the closed ball B(obj,range_radius)
        """
        result = []
        queue = deque([self.root])
        while queue:
            node = queue.pop()
            if rect.contains_rect(node.rect):
                result.extend(get_points(node))

            if node.is_leaf:
                result.extend((p for p in node.entry_list if rect.contains(p)))
            else:
                queue.extend([n for n in node.subtree if rect.does_intersect(n.rect)])
        return result

    def average_leaf_volume(self):
        sum = 0
        counter = 0
        queue = deque([self.root])
        while queue:
            node = queue.pop()
            if node.is_leaf:
                sum += node.rect.get_volume()
                counter += 1
            else:
                queue.extend(node.subtree)
        return sum/counter, counter


    def get_volume(self):
        sum = 0
        queue = deque([self.root])
        while queue:
            node = queue.pop()
            sum += node.rect.get_volume()
            if not node.is_leaf:
                queue.extend(node.subtree)
        return sum


    def get_all_points(self):
        return get_all_points(self.root)

    def get_depth(self):
        return self.root.get_depth()

    def get_max_inner_node_branching(self):
        return self.root.get_max_inner_node_branching()

class InnerNode:

    def __init__(self, tree_ptr, rect, subtree):
        self.tree_ptr = tree_ptr
        self.rect = rect
        self.subtree = subtree
        self.is_leaf = False
        self.volume = rect.get_volume()


    def insert(self, obj):
        insert_node = self.choose_subtree(obj, self.subtree)
        insert_node.rect = get_enlarged_rect(insert_node.rect, obj)
        split_node = insert_node.insert(obj)


        if split_node != None:
            # if not self.rect.contains_rect(insert_node.rect):
            #     print("Self not containing insert node rect")
            # if not self.rect.contains_rect(split_node.rect):
            #     print("Self not containing split node rect")
            if self.subtree.__len__() < self.tree_ptr.max_inner_node_children:
                self.subtree.append(split_node)
                return None
            else:
                return self.split(split_node)


    def choose_subtree(self, obj, subtree_list):
        '''
        if all([n.is_leaf for n in subtree_list]):
            rect_list = [n.rect for n in subtree_list]
            overlap = get_overlap(rect_list)
            minimal_enlargement_of_overlap = math.inf
            minimal_area_enlargement = math.inf
            choose_node = None
            for ind, node in enumerate(subtree_list):
                enl_rect = get_enlarged_rect(node.rect, obj)
                #Create new lists:
                new_rect_list = rect_list.copy()
                new_rect_list[ind] = enl_rect
                #Calculate enlargement of overlap:
                enl_of_overlap = get_overlap_single(new_rect_list, ind)
                if enl_of_overlap - overlap[ind] < minimal_enlargement_of_overlap:
                    choose_node = node
                    minimal_enlargement_of_overlap = enl_of_overlap - overlap[ind]
            return choose_node
        else:
        '''
        minimal_area_enlargement = math.inf
        choose_node = None
        for ind, node in enumerate(subtree_list):
            new_vol = get_enlarged_rect(node.rect, obj).get_volume()
            if new_vol - node.volume < minimal_area_enlargement:
                minimal_area_enlargement = new_vol - node.volume
                choose_node = node
        return choose_node


    def split(self, node):
        self.subtree.append(node)
        node_1, node_2 = self.pick_seeds()
        self.subtree.remove(node_1)
        self.subtree.remove(node_2)
        group_1, group_2 = [node_1], [node_2]
        bounding_rect_group_1 = node_1.rect
        bounding_rect_group_2 = node_2.rect
        max_size = self.tree_ptr.max_inner_node_children - self.tree_ptr.min_inner_node_children+1
     #   print("is correct before", [n.is_correct() for n in self.subtree])
        while self.subtree != [] and group_1.__len__() < max_size and group_2.__len__() < max_size:
            next_node = RTree.get_next_node(self.subtree, bounding_rect_group_1, bounding_rect_group_2)
            #Find group to insert node
            k = RTree.get_group_to_insert_node(next_node, bounding_rect_group_1, bounding_rect_group_2)
            if k == 1:
                group_1.append(next_node)
                bounding_rect_group_1 = get_union(bounding_rect_group_1, next_node.rect)
            else:
                group_2.append(next_node)
                bounding_rect_group_2 = get_union(bounding_rect_group_2, next_node.rect)
            self.subtree.remove(next_node)

        if self.subtree != []:
            if group_1.__len__() == max_size:
                group_2.extend(self.subtree)
                rect_list = [bounding_rect_group_2]
                rect_list.extend([n.rect for n in self.subtree])
                bounding_rect_group_2 = get_union_list(rect_list)
            else:
                group_1.extend(self.subtree)
                rect_list = [bounding_rect_group_1]
                rect_list.extend([n.rect for n in self.subtree])
                bounding_rect_group_1 = get_union_list(rect_list)
        #Update current node:
        self.rect = bounding_rect_group_1
        self.subtree = group_1
        #Return the second node for insertion
        return InnerNode(self.tree_ptr, bounding_rect_group_2, group_2)


    def pick_seeds(self):
        node_1, node_2 = None, None
        greatest_common_rect_area = -math.inf
        n = self.subtree.__len__()
        for i in range(n):
            for j in range(i):
                n_1, n_2 = self.subtree[i], self.subtree[j]
                vol = get_union(n_1.rect, n_2.rect).get_volume()
                d = vol - n_1.volume - n_2.volume
                if d > greatest_common_rect_area:
                    node_1, node_2 = n_1, n_2
                    greatest_common_rect_area = d
        return node_1, node_2


    def is_correct(self):
        for n in self.subtree:
            if not n.is_correct():
                return False
            if not self.rect.contains_rect(n.rect):
                #print("Inner Node not correct", self.rect.max_dimensions, n.rect.max_dimensions)
                return False
        return True


    def get_depth(self):
        return max([n.get_depth() for n in self.subtree])+1

    def get_max_inner_node_branching(self):
        r = max([n.get_max_inner_node_branching() for n in self.subtree])
        return max(r, self.subtree.__len__())

class LeafNode:

    def __init__(self, tree_ptr, rect, entry_list):
        self.tree_ptr = tree_ptr
        self.rect = rect
        self.entry_list = entry_list
        self.is_leaf = True
        self.volume = rect.get_volume()


    def insert(self, obj):
        if self.entry_list.__len__() < self.tree_ptr.max_leaf_node_elements:
            self.entry_list.append(obj)
            return None
        else:
            return self.split(obj)


    def split(self, obj):
        self.entry_list.append(obj)
        obj_1, obj_2 = self.pick_seeds()
        group_1, group_2 = [obj_1], [obj_2]
        bounding_rect_group_1 = Rect(obj_1, obj_1)
        bounding_rect_group_2 = Rect(obj_2, obj_2)
        max_size = self.tree_ptr.max_leaf_node_elements - self.tree_ptr.min_leaf_node_elements+1
        while self.entry_list != [] and group_1.__len__() < max_size and group_2.__len__() < max_size:
            next_obj = RTree.get_next_obj(self.entry_list, bounding_rect_group_1, bounding_rect_group_2)
            #Find group to insert node
            k = RTree.get_group_to_insert_obj(next_obj, bounding_rect_group_1, bounding_rect_group_2)
            if k == 1:
                group_1.append(next_obj)
                bounding_rect_group_1 = get_enlarged_rect(bounding_rect_group_1, next_obj)
            else:
                group_2.append(next_obj)
                bounding_rect_group_2 = get_enlarged_rect(bounding_rect_group_2, next_obj)
            self.entry_list.remove(next_obj)

        if self.entry_list != []:
            if group_1.__len__() == max_size:
                group_2.extend(self.entry_list)
                bounding_rect_group_2 = get_union(bounding_rect_group_2, get_union_points_list(self.entry_list))
            else:
                group_1.extend(self.entry_list)
                bounding_rect_group_1 = get_union(bounding_rect_group_1, get_union_points_list(self.entry_list))
        #Update current node:
        self.rect = bounding_rect_group_1
        self.entry_list = group_1
        #Return the second node for insertion
        return LeafNode(self.tree_ptr, bounding_rect_group_2, group_2)


    def pick_seeds(self):
        obj_1, obj_2 = None, None
        greatest_common_rect_area = -math.inf
        n = self.entry_list.__len__()
        for i in range(n):
            for j in range(i):
                e_1, e_2 = self.entry_list[i], self.entry_list[j]
                vol = get_union_points(e_1, e_2).get_volume()
                d = vol
                if d > greatest_common_rect_area:
                    obj_1, obj_2 = e_1, e_2
                    greatest_common_rect_area = d
        return obj_1, obj_2


    def is_correct(self):
        for p in self.entry_list:
            if not self.rect.contains(p):
                #print("Leaf Node not correct")
                return False
        return True

    def get_depth(self):
        return 1

    def get_max_inner_node_branching(self):
        return 0


def bulk_load_r_tree_morton(point_list, max_number_of_leaf_points, min_number_of_leaf_points, max_number_of_inner_nodes, min_number_of_inner_nodes, sq_metric):
    n = point_list.__len__()
    m = max_number_of_leaf_points
    l = max_number_of_inner_nodes
    r_tree = RTree(None, max_number_of_inner_nodes, min_number_of_inner_nodes, max_number_of_leaf_points, min_number_of_leaf_points, sq_metric)
    point_list.sort(key=lambda x:morton_code(x))
    node_list = []
    for i in range(0,n,m):
        current_list = point_list[i:i+m]
        rect = get_union_points_list(current_list)
        node_list.append(LeafNode(r_tree, rect, current_list))
    while node_list.__len__() > 1:
        node_list.sort(key=lambda x:morton_code(x.rect.get_center()))
        new_node_list = []
        w = node_list.__len__()
        for i in range(0, w, l):
            current_list = node_list[i:i + l]
            rect = get_union_list([nd.rect for nd in current_list])
            new_node_list.append(InnerNode(r_tree, rect, current_list))
        node_list = new_node_list
    r_tree.root = node_list[0]
    return r_tree


def bulk_load_r_tree_str_2(point_list, max_number_of_leaf_points, min_number_of_leaf_points, max_number_of_inner_nodes, min_number_of_inner_nodes):
    r = point_list.__len__()
    n = max_number_of_leaf_points
    l = max_number_of_inner_nodes
    r_tree = RTree(None, max_number_of_inner_nodes, min_number_of_inner_nodes, max_number_of_leaf_points, min_number_of_leaf_points)
    point_list.sort(key=lambda x:x[0])
    p = math.ceil(r/n)
    s = math.ceil(math.sqrt(p))
    node_list = []
    for i in range(0,r,s*n):
        current_list = point_list[i:i+s*n]
        current_list.sort(key=lambda x:x[1])
        for j in range(0,current_list.__len__(),n):
            rect = get_union_points_list(current_list[j:j+n])
            node_list.append(LeafNode(r_tree, rect, current_list[j:j+n]))
    while node_list.__len__() > 1:
        r = node_list.__len__()
        p = math.ceil(r / l)
        s = math.ceil(math.sqrt(p))
        node_list.sort(key=lambda x: x.rect.get_center()[0])
        new_node_list = []
        for i in range(0, r, s*l):
            current_list = node_list[i:i + s*l]
            current_list.sort(key=lambda x:x.rect.get_center()[1])
            for j in range(0,current_list.__len__(),l):
                rect = get_union_list([nd.rect for nd in current_list[j:j+l]])
                new_node_list.append(InnerNode(r_tree, rect, current_list[j:j+l]))
        node_list = new_node_list
    r_tree.root = node_list[0]
    return r_tree


def bulk_load_r_tree_str(point_list, max_number_of_leaf_points, min_number_of_leaf_points, max_number_of_inner_nodes, min_number_of_inner_nodes, sq_metric):
    r_tree = RTree(None, max_number_of_inner_nodes, min_number_of_inner_nodes, max_number_of_leaf_points, min_number_of_leaf_points, sq_metric)
    dim = point_list[0].__len__()
    node_list = slice_data_points(point_list, 0, dim, dim, max_number_of_leaf_points, r_tree)
    while node_list.__len__() > 1:
        node_list = slice_data(node_list, 0, dim, dim, max_number_of_inner_nodes, r_tree)
    r_tree.root = node_list[0]
    return r_tree



def slice_data(node_list, coord, dim, current_dim, max_number_of_children, r_tree):
    node_list.sort(key=lambda x: x.rect.get_center()[coord])
    if current_dim == 1:
        # This is the base case
        result = []
        for i in range(0,node_list.__len__(),max_number_of_children):
            rect = get_union_list([nd.rect for nd in node_list[i:i+max_number_of_children]])
            result.append(InnerNode(r_tree, rect, node_list[i:i+max_number_of_children]))
        return result
    else:
        #Divide the node list in slabs:
        number_nodes_in_slabs = max_number_of_children*math.ceil(math.pow(math.ceil(node_list.__len__()/max_number_of_children),(current_dim-1)/current_dim))
        result = []
        for i in range(0,node_list.__len__(), number_nodes_in_slabs):
            slab = node_list[i:i+number_nodes_in_slabs]
            #Recursively slice the slabs
            result.extend(slice_data(slab, coord+1, dim, current_dim - 1, max_number_of_children, r_tree))
        return result


def slice_data_points(point_list, coord, dim, current_dim, max_number_of_children, r_tree):
    point_list.sort(key=lambda x: x[coord])
    if current_dim == 1:
        # This is the base case
        result = []
        for i in range(0,point_list.__len__(),max_number_of_children):
            rect = get_union_points_list(point_list[i:i+max_number_of_children])
            result.append(LeafNode(r_tree, rect, point_list[i:i+max_number_of_children]))
        return result
    else:
        #Divide the node list in slabs:
        number_nodes_in_slabs = max_number_of_children*math.ceil(math.pow(math.ceil(point_list.__len__()/max_number_of_children),(current_dim-1)/current_dim))
        result = []
        for i in range(0,point_list.__len__(), number_nodes_in_slabs):
            slab = point_list[i:i+number_nodes_in_slabs]
            #Recursively slice the slabs
            result.extend(slice_data_points(slab, coord+1, dim, current_dim - 1, max_number_of_children, r_tree))
        return result

