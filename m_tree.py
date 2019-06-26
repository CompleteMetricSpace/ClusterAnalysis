import math, time, random
from collections import deque
from morton import morton_code



node_point_dict = {}
def get_points(node):
    """
    :param node: a MTree-Node
    :return: The list of all points in this node
    """
    if node in node_point_dict:
        return node_point_dict[node]
    points = []
    queue = deque([node])
    while queue:
        n = queue.pop()
        if n.is_leaf:
            points.extend([p[0] for p in n.entry_list])
        else:
            queue.extend(n.subtree)
    node_point_dict[node] = points
    return points


class MTree:
    """
    M-Tree implementation as described in
    "P. Zezula, P. Ciaccia, and F. Rabitti. M-tree: A dynamic index for similarity queries in multimedia databases",
    Technical Report 7, HERMES ESPRIT LTR Project, 1996

    This implementation uses a modified m_RAD splitting policy

    """
    def __init__(self, metric, max_number_of_children=50, max_number_of_points_in_leafs = 100, radius_error=1e-7, sq_metric = None):
        self.root = None
        self.max_number_of_children = max_number_of_children
        self.max_number_of_points_in_leafs = max_number_of_points_in_leafs
        self.metric = metric
        self.radius_error = radius_error
        if sq_metric == None:
            self.sq_metric = lambda x,y: metric(x,y)**2
        else:
            self.sq_metric = sq_metric

    def insert(self, obj):
        """
        Inserts an object in this M-Tree
        :param obj: an object
        :return: Nothing is returned
        """
        if self.root is None:
            #If root is not initialized, then initialize it as a leaf
            node = LeafNode(self, obj, 0, None, [(obj, 0)])
            self.root = node
        else:
            #Else insert the object into the root node and increase its radius
            r = self.metric(self.root.obj, obj)
            if self.root.covering_radius < r:
                self.root.covering_radius = r + self.radius_error
            split_node = self.root.insert(obj)
            if split_node != None:
                #If the root node splits, then initalize a new node (that is the future root node) with
                #routing object as the current object
                split_node_1 = split_node[0]
                split_node_2 = split_node[1]
                parent = self.root.obj
                #Set split node distance to (future) parent node
                split_node_2.distance_to_parent = self.metric(parent, split_node_2.obj)
                #Set new covering radius as the maximum
                r = max(self.root.covering_radius, split_node_2.distance_to_parent + split_node_2.covering_radius)
                split_node_1.distance_to_parent = 0
                #Create a new inner node and promote it to be the root
                self.root = InnerNode(self, parent, r, None, [split_node_1, split_node_2])


    def get_points_in_range(self, obj, range_radius):
        result = []
        d = self.metric(obj, self.root.obj)
        self.root.get_points_in_range(obj, range_radius, d, result)
        return result


    def get_points_in_range_it(self, obj, range_radius):
        """
        Get all points that are included in the closed ball of radius range_radius around obj
        This is the iterative method
        :param obj: a point in space (query point)
        :param range_radius: a number > 0 (query radius)
        :return: the list of all points that are in B(obj, range_radius)
        """
        # We never use the metric, since we assume that sq_metric is always faster
        result = []
        dist = self.sq_metric(self.root.obj, obj)
        queue = deque([(self.root, dist)])
        while queue:
            node, distance_to_this_node_sq = queue.pop()
            if node.is_leaf:
                for (p, r) in node.entry_list:
                    #This checks if distance_to_this_node + r <= range+radius
                    #By the triangle inequality we can then safely add this point
                    if range_radius >= r and distance_to_this_node_sq <= (range_radius-r)*(range_radius-r):
                        result.append(p)
                        continue
                    #This checks if abs(distance_to_this_node - r) <= range_radius, since if this
                    #is not the case, then we can safely discard the point p
                    d = distance_to_this_node_sq
                    r_1 = range_radius+r
                    r_2 = r-range_radius
                    if (d >= r*r and d <= r_1*r_1) or (d < r*r and (r_2 <= 0 or d >= r_2*r_2)):
                        #This checks whether the point p is really in the range_radius of obj
                        if self.sq_metric(obj, p) <= range_radius*range_radius:
                            result.append(p)
            else:
                for subnode in node.subtree:
                    sum_of_radii = range_radius + subnode.covering_radius
                    #This checks if abs(distance_to_this_node - subnode.distance_to_parent) <= sum_of_radii, since if this
                    #is not the case, then we can safely discard the subnode
                    dn = distance_to_this_node_sq
                    dp = subnode.distance_to_parent
                    dp_sq = subnode.distance_to_parent*subnode.distance_to_parent
                    srdp_sq = (sum_of_radii+subnode.distance_to_parent)*(sum_of_radii+subnode.distance_to_parent)
                    dpsr_sq = (subnode.distance_to_parent-sum_of_radii)*(subnode.distance_to_parent-sum_of_radii)
                    if (dn >= dp_sq and dn <= srdp_sq) or (dn<=dp_sq and (dp <=sum_of_radii or dpsr_sq <= dn)):
                        d = self.sq_metric(subnode.obj, obj)
                        #Check if distance to this subnode is small enough. If this is not
                        #the case, then we can savely discard this node
                        if d <= sum_of_radii*sum_of_radii:
                            #Check if this subnode is completely included in the range ball around obj
                            #If this is true, then we can just add all points of this subnode
                            if range_radius >= subnode.covering_radius and \
                                d <= (range_radius - subnode.covering_radius)*(range_radius - subnode.covering_radius):
                                result.extend(get_points(subnode))
                            else:
                                queue.append((subnode, d))
        return result


    def get_points_in_range_new(self, range_radius):
        result = {x:[] for x in self.get_all_points()}
        leafes = self.get_leaf_nodes()
        leafes_points = {n:[p[0] for p in n.entry_list] for n in leafes}
        range_sq = range_radius*range_radius
        n = leafes.__len__()
        for i in range(n):
            print(i)
            for j in range(i+1):
                p = leafes[i]
                q = leafes[j]
                c = p.covering_radius + q.covering_radius+range_radius
                d = self.sq_metric(p.obj, q.obj)
                if range_radius >= p.covering_radius + q.covering_radius and d <= (range_radius-p.covering_radius-q.covering_radius)**2:
                    for obj in leafes_points[p]:
                        result[obj].extend(leafes_points[q])
                    for obj in leafes_points[q]:
                        result[obj].extend(leafes_points[p])
                    continue
                if d <= c:
                    for obj_1, r_1 in p.entry_list:
                        d2 = self.sq_metric(obj_1, q.obj)
                        if q.covering_radius <= range_radius and d2 <= (range_radius-q.covering_radius)**2:
                            result[obj_1].extend(leafes_points[q])
                            for obj_2, r_2 in q.entry_list:
                                result[obj_2].append(obj_1)
                            continue
                        for obj_2, r_2 in q.entry_list:
                            if self.sq_metric(obj_1, obj_2) <=range_sq:
                                result[obj_1].append(obj_2)
                                result[obj_2].append(obj_1)
        return result


    def get_all_points(self):
        result = []
        self.root.append_all_points(result)
        return result


    def get_leaf_nodes(self):
        result = []
        queue = deque([self.root])
        while queue:
            node = queue.pop()
            if node.is_leaf:
                result.append(node)
            else:
                queue.extend(node.subtree)
        return result

    def get_volume(self, dim):
        return self.root.get_volume(dim)

    def is_correct(self):
        return self.root.is_correct()

class InnerNode:

    def __init__(self, tree_ptr, obj, covering_radius, distance_to_parent, subtree):
        self.tree_ptr = tree_ptr
        self.obj = obj
        self.covering_radius = covering_radius
        self.distance_to_parent = distance_to_parent
        self.subtree = subtree
        self.number_of_points = sum([n.number_of_points for n in subtree])
        self.is_leaf = False

    def split(self, node):
        self.subtree.append(node)
        # Determine the objects that should be promoted: (this is the modified m_RAD - strategy)
        # We try all pairs of nodes and determine which pair does inflict the least
        # maximum of radii (sum of radii in the original paper)
        parent_1 = None #Is a node
        parent_1_subtree_list = []
        parent_1_covering_radius = 0
        parent_2 = None #Is a node
        parent_2_subtree_list = []
        parent_2_covering_radius = 0
        sum_of_radii = math.inf
        length = self.subtree.__len__()
        for i in range(length):
            for j in range(i):
                node_1 = self.subtree[i]
                node_2 = self.subtree[j]
                if node_1 is node_2:
                    continue
                # Do generalized hyperplane split
                subtree_list_1 = []
                subtree_list_2 = []
                covering_radius_1 = 0
                covering_radius_2 = 0

                # Redistribute the nodes in subtree among the two new nodes
                for node in self.subtree:
                    distance_to_node_1 = self.tree_ptr.metric(node_1.obj, node.obj)
                    distance_to_node_2 = self.tree_ptr.metric(node_2.obj, node.obj)
                    if distance_to_node_2 < distance_to_node_1:
                        subtree_list_2.append(node)
                        # Set node distance to parent and the new covering radius
                        covering_radius_2 = max(covering_radius_2, distance_to_node_2 + node.covering_radius)
                    else:
                        subtree_list_1.append(node)
                        covering_radius_1 = max(covering_radius_1, distance_to_node_1 + node.covering_radius)

                if max(covering_radius_1, covering_radius_2) < sum_of_radii:
                    #Save node-pair!
                    parent_1 = node_1
                    parent_2 = node_2
                    parent_1_subtree_list = subtree_list_1
                    parent_2_subtree_list = subtree_list_2
                    parent_1_covering_radius = covering_radius_1
                    parent_2_covering_radius = covering_radius_2
                    sum_of_radii = max(covering_radius_1,  covering_radius_2)

        #Set up the current node
        self.obj = parent_1.obj
        self.covering_radius = parent_1_covering_radius + self.tree_ptr.radius_error
        self.subtree = parent_1_subtree_list
        # Set distances:
        for node in self.subtree:
            node.distance_to_parent = self.tree_ptr.metric(node.obj, parent_1.obj)
        self.number_of_points = sum([x.number_of_points for x in self.subtree])

        #Set up second node
        for node in parent_2_subtree_list:
            node.distance_to_parent = self.tree_ptr.metric(node.obj, parent_2.obj)

        # Distance to parent is set later
        return self, InnerNode(self.tree_ptr, parent_2.obj, parent_2_covering_radius + self.tree_ptr.radius_error, None, parent_2_subtree_list)


    def insert(self, obj):
        """
        Inserts a new object in this node
        :param obj: an object
        :return: a node if this node splits and nothing otherwise
        """
        self.number_of_points += 1
        node_not_enlarging_covering_radius = None
        distance_to_node_not_enlarging_cr = math.inf
        node_with_least_increase_of_covering_radius = None
        distance_to_node_with_least_increase_of_cr = math.inf
        least_increase_of_cr = math.inf
        #First look at the nodes that do not require an increase of covering radius
        #If such nodes exist, then pick the one that has a routing object that is most close to obj
        #If all nodes require an increase of covering radius, pick the one requiring the least increase
        for node in self.subtree:
            r = self.tree_ptr.metric(node.obj, obj)
            if r <= node.covering_radius:
                if distance_to_node_not_enlarging_cr > r:
                    node_not_enlarging_covering_radius = node
                    distance_to_node_not_enlarging_cr = r
            else:
                if node_not_enlarging_covering_radius == None:
                    if r - node.covering_radius < least_increase_of_cr:
                        node_with_least_increase_of_covering_radius = node
                        least_increase_of_cr = r - node.covering_radius
                        distance_to_node_with_least_increase_of_cr = r
        if node_not_enlarging_covering_radius != None:
            insertion_node = node_not_enlarging_covering_radius
        else:
            insertion_node = node_with_least_increase_of_covering_radius
            insertion_node.covering_radius = distance_to_node_with_least_increase_of_cr + self.tree_ptr.radius_error



        #Insert the object into the node
        split_node = insertion_node.insert(obj)

        if split_node == None:
            return None
        else:
            split_node_1 = split_node[0]
            split_node_2 = split_node[1]
            #Set the distance of the split node, since this is not done right away in the split method
            split_node_1.distance_to_parent = self.tree_ptr.metric(self.obj, split_node_1.obj)
            split_node_2.distance_to_parent = self.tree_ptr.metric(self.obj, split_node_2.obj)
            # TODO: Should the covering radius be increased here?
            # No, since the split node still is completely included in this node
            if self.subtree.__len__() < self.tree_ptr.max_number_of_children:
                #Just append the splitting node 2, since the split node 1 is already appended
                self.subtree.append(split_node_2)
                return None
            else:
                # Split this node in this case
                #Just pass it further up
                return self.split(split_node_2)


    def get_points_in_range(self, obj, range_radius, distance_to_this_node, result):
        for node in self.subtree:
            sum_of_radii = range_radius + node.covering_radius
            if abs(distance_to_this_node-node.distance_to_parent) <= sum_of_radii+self.tree_ptr.radius_error:
                d = self.tree_ptr.metric(node.obj, obj)
                if d<= sum_of_radii:
                    if d + node.covering_radius <= range_radius:
                        node.append_all_points(result)
                    else:
                        node.get_points_in_range(obj, range_radius, d, result)


    def get_number_of_points_in_range(self, obj: tuple, range_radius: float, distance_to_this_node: float) -> float:
        sum = 0
        for node in self.subtree:
            sum_of_radii = range_radius + node.covering_radius
            if abs(distance_to_this_node-node.distance_to_parent) <= sum_of_radii+self.tree_ptr.radius_error:
                d = self.tree_ptr.metric(node.obj, obj)
                if d<= sum_of_radii:
                    if d + node.covering_radius <= range_radius:
                        sum += node.number_of_points
                    else:
                        sum += node.get_number_of_points_in_range(obj, range_radius, d)
        return sum


    def append_all_points_set(self, point_set, result_map):
        for node in self.subtree:
            node.append_all_points_set(point_set, result_map)

    def append_all_points(self, point_list):
        for node in self.subtree:
            node.append_all_points(point_list)


    def get_points(self):
        res = []
        for n in self.subtree:
            res.extend(n.get_points())
        return res


    def get_volume(self, dim):
        volume = math.pow(math.pi, dim/2)/math.gamma(dim/2+1)*math.pow(self.covering_radius, dim)
        for x in self.subtree:
            volume += x.get_volume(dim)
        return volume


    def is_correct(self):
        for n in self.subtree:
            if not n.is_correct():
                return False
            if n.distance_to_parent > self.covering_radius+self.tree_ptr.radius_error:
                return False
        return True



class LeafNode:

    def __init__(self, tree_ptr, obj, covering_radius, distance_to_parent, entry_list):
        self.tree_ptr = tree_ptr
        self.obj = obj
        self.covering_radius = covering_radius
        self.distance_to_parent = distance_to_parent
        self.entry_list = entry_list
        self.number_of_points = entry_list.__len__()
        self.is_leaf = True

    def insert(self, obj):

        self.number_of_points += 1
        if self.entry_list.__len__() >= self.tree_ptr.max_number_of_points_in_leafs:
            #Just split
            return self.split(obj)
        else:
            #Just append and optionally increase the covering radius
            r = self.tree_ptr.metric(self.obj, obj)
            if self.covering_radius < r:
                self.covering_radius = r + self.tree_ptr.radius_error
            self.entry_list.append((obj, r))
            return None


    def split(self, obj):
        r = self.tree_ptr.metric(self.obj, obj)
        self.entry_list.append((obj, r))
        # Determine the farthest point from self.obj:
        # parent_1 = self.obj
        parent_2, r_max = max(self.entry_list, key=lambda x: x[1])
        # Do generalized hyperplane split
        entry_list_1 = []
        entry_list_2 = []
        covering_radius_1 = 0
        covering_radius_2 = 0
        for (elem, dist) in self.entry_list:
            r = self.tree_ptr.metric(elem, parent_2)
            if r < dist:
                entry_list_2.append((elem, r))
                covering_radius_2 = max(covering_radius_2, r)
            else:
                entry_list_1.append((elem, dist))
                covering_radius_1 = max(covering_radius_1, dist)

        self.entry_list = entry_list_1
        self.covering_radius = covering_radius_1 + self.tree_ptr.radius_error
        self.number_of_points = entry_list_1.__len__()

        # Distance to parent is set later
        return self, LeafNode(self.tree_ptr, parent_2, covering_radius_2 + self.tree_ptr.radius_error, None, entry_list_2)


    def get_points_in_range(self, obj, range_radius, distance_to_this_node, result):
        for (p, r) in self.entry_list:
            if abs(distance_to_this_node - r) <= range_radius+self.tree_ptr.radius_error:
                if self.tree_ptr.metric(obj, p) <= range_radius:
                    result.append(p)


    def get_number_of_points_in_range(self, obj: tuple, range_radius: float, distance_to_this_node: float) -> int:
        sum = 0
        for (p, r) in self.entry_list:
            if abs(distance_to_this_node - r) <= range_radius+self.tree_ptr.radius_error:
                if self.tree_ptr.metric(obj, p) <= range_radius:
                    sum += 1
        return sum


    def append_all_points_set(self, point_set, result_map):
        for p in point_set:
            result_map[p].extend([q[0] for q in self.entry_list])

    def append_all_points(self, point_list):
        for (p,r) in self.entry_list:
            point_list.append(p)


    def get_points(self):
        return list(map(lambda x:x[0], self.entry_list))

    def get_volume(self, dim):
        return math.pow(math.pi, dim/2)/math.gamma(dim/2+1)*math.pow(self.covering_radius, dim)

    def is_correct(self):
        for (p,r) in self.entry_list:
            if r > self.covering_radius+self.tree_ptr.radius_error:
                return False
        return True





def get_union_points_list(point_list, sq_metric):
    n = point_list.__len__()
    average_point = tuple(map(lambda x:sum(x)/n, zip(*point_list)))
    middle_point = None
    distance_to_average_point = math.inf
    for p in point_list:
        d = sq_metric(p, average_point)
        if d < distance_to_average_point:
            middle_point = p
            distance_to_average_point = d
    radius = 0
    child_list = []
    for p in point_list:
        d = math.sqrt(sq_metric(p, middle_point))
        if d > radius:
            radius = d
        child_list.append((p, d))
    return middle_point, radius, child_list


def get_union_list(ball_list, sq_metric):
    n = ball_list.__len__()
    average_point = tuple(map(lambda x:sum(x)/n, zip(*tuple(b[0] for b in ball_list))))
    middle_point = None
    distance_to_average_point = math.inf
    for p, r in ball_list:
        d = sq_metric(p, average_point)
        if d < distance_to_average_point:
            middle_point = p
            distance_to_average_point = d
    radius = 0
    child_dist_list = []
    for p, r in ball_list:
        d = math.sqrt(sq_metric(p, middle_point))
        if d+r > radius:
            radius = d+r
        child_dist_list.append(d)
    return middle_point, radius, child_dist_list



def bulk_load_m_tree_morton(point_list, max_number_of_leaf_points, max_number_of_inner_nodes, metric, sq_metric):
    n = point_list.__len__()
    m = max_number_of_leaf_points
    l = max_number_of_inner_nodes
    m_tree = MTree(metric,max_number_of_children=max_number_of_inner_nodes, max_number_of_points_in_leafs=max_number_of_leaf_points, sq_metric=sq_metric)
    point_list.sort(key=lambda x:morton_code(x))
    node_list = []
    for i in range(0,n,m):
        current_list = point_list[i:i+m]
        obj, radius, children = get_union_points_list(current_list, sq_metric)
        node_list.append(LeafNode(m_tree, obj, radius, None, children))
    while node_list.__len__() > 1:
        node_list.sort(key=lambda x:morton_code(x.obj))
        new_node_list = []
        w = node_list.__len__()
        for i in range(0, w, l):
            current_list = node_list[i:i + l]
            obj, radius, children_distances = get_union_list([(nd.obj, nd.covering_radius) for nd in current_list], sq_metric)
            for j in range(current_list.__len__()):
                current_list[j].distance_to_parent = children_distances[j]
            new_node_list.append(InnerNode(m_tree, obj, radius, None, current_list))
        node_list = new_node_list
    m_tree.root = node_list[0]
    return m_tree


def bulk_load_m_tree_str_2(point_list, max_number_of_leaf_points, max_number_of_inner_nodes, metric, sq_metric):
    r = point_list.__len__()
    n = max_number_of_leaf_points
    l = max_number_of_inner_nodes
    m_tree = MTree(metric, max_number_of_inner_nodes, max_number_of_leaf_points, sq_metric=sq_metric)

    point_list.sort(key=lambda x:x[0])
    p = math.ceil(r/n)
    s = math.ceil(math.sqrt(p))
    node_list = []
    for i in range(0,r,s*n):
        current_list = point_list[i:i+s*n]
        current_list.sort(key=lambda x:x[1])
        for j in range(0,current_list.__len__(),n):
            obj, radius, children = get_union_points_list(current_list[j:j+n], sq_metric)
            node_list.append(LeafNode(m_tree, obj, radius, None, children))
    while node_list.__len__() > 1:
        r = node_list.__len__()
        p = math.ceil(r / l)
        s = math.ceil(math.sqrt(p))
        node_list.sort(key=lambda x: morton_code(x.obj))
        new_node_list = []
        for i in range(0, r, s*l):
            current_list = node_list[i:i + s*l]
            current_list.sort(key=lambda x:x.obj[1])
            for j in range(0,current_list.__len__(),l):
                current_list_this = current_list[j:j+l]
                obj, radius, children_distances = get_union_list([(nd.obj, nd.covering_radius) for nd in current_list_this],sq_metric)
                new_node_list.append(InnerNode(m_tree, obj, radius, None, current_list_this))
                for k in range(current_list_this.__len__()):
                    current_list_this[k].distance_to_parent = children_distances[k]
        node_list = new_node_list
    m_tree.root = node_list[0]
    return m_tree
