import csv
import math
import random
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import colorsys
from collections import deque
import time
#from MTree import MTree, LeafNode, InnerNode
import m_tree

from m_tree import MTree, LeafNode, InnerNode, bulk_load_m_tree_morton, bulk_load_m_tree_str_2
from kd_tree import KDTree, create_kd_tree
import kd_tree as kdt
import cProfile as profile
import metrics as mtr
from list_of_clusters import ClusterList, ClusterNode, build_cluster_list_it
from morton import morton_code
from r_tree import RTree, bulk_load_r_tree_morton, bulk_load_r_tree_str, Rect

name = '/home/andrei/Documents/University/Computerpraktikum/Python/ExtractedFiles/Cluster/cod-rna.5000.csv'
point_list = []
with open(name) as csvfile:
    data = csv.reader(csvfile)
    for row in data:
        point_list.append(tuple([float(s) for s in row]))
ref_point = point_list[0]
point_list = point_list[0:10000]
dimension = ref_point.__len__()
use_sq_metric = mtr.choose_sq_metric(dimension)
#point_list.sort(key=lambda x: morton_code(x))


def test_m_tree(tau, dim=2, sample_size=2000, sq_metric=mtr.sq_dist_2, profiler=False, leaf_points=20, inner_children=6):
    print("---- Testing M-Tree range search for %s sample points" % sample_size, "and tau = ", tau, " --- ")
    start_time = time.time()
    m_tree = MTree(mtr.dist_mem, inner_children, leaf_points, sq_metric=sq_metric)
    for p in point_list:
        m_tree.insert(p)
    print("--- %s seconds " % (time.time() - start_time), " for building M-Tree ---")
    print("M-Tree is correct: ", m_tree.is_correct())
    print("M-Tree has volume:", m_tree.get_volume(dim))
    leafs = m_tree.get_leaf_nodes()
    print("M-Tree has number of leaf-nodes: ", leafs.__len__())
    print("M-Tree leaf node has average covering radius of", sum([l.covering_radius for l in leafs])/leafs.__len__())
    result = {}
    smpl = random.sample(point_list, sample_size)
    start_time = time.time()
    pr = profile.Profile()
    if profiler:
        pr.enable()
    for p in smpl:
        result[p] = m_tree.get_points_in_range_it(p, tau)
    if profiler:
        pr.disable()
    print("--- %s seconds " % ((time.time() - start_time)/sample_size), " per point query")
    if profiler:
        pr.dump_stats('m_tree_min_rad.pstat')
    return result


def test_kd_tree(tau, dim=2, sample_size=2000, sq_metric=mtr.sq_dist_2, profiler=False, leaf_points=100):
    print("---- Testing KD-Tree range search for %s sample points" % sample_size, " and tau = ", tau, " --- ")
    start_time = time.time()
    kd_tree = create_kd_tree(point_list, leaf_points, sq_metric)
    print("--- %s seconds " % (time.time() - start_time), " for building KD-Tree ---")
    result = {}
    smpl = random.sample(point_list, sample_size)
    start_time = time.time()
    pr = profile.Profile()
    if profiler:
        pr.enable()
    for p in smpl:
        result[p] = kd_tree.get_points_in_range_it(p, tau)
    if profiler:
        pr.disable()
    print("--- %s seconds " % ((time.time() - start_time)/sample_size), " per point query")
    if profiler:
        pr.dump_stats('kd_tree.pstat')
    return result


def test_cluster_list(tau, dim=2, sample_size=2000, sq_metric=mtr.sq_dist_2, profiler=False, disk_radius=0.075):
    print("---- Testing Cluster-List range search for %s sample points" % sample_size, " and tau = ", tau, " --- ")
    start_time = time.time()
    cl = build_cluster_list_it(point_list, disk_radius , sq_metric)
    print("--- %s seconds " % (time.time() - start_time), " for building Cluster-List ---")
    result = {}
    smpl = random.sample(point_list, sample_size)
    start_time = time.time()
    pr = profile.Profile()
    if profiler:
        pr.enable()
    for p in smpl:
        result[p] = cl.search(p, tau)
    if profiler:
        pr.disable()
    print("--- %s seconds " % ((time.time() - start_time)/sample_size), " per point query")
    if profiler:
        pr.dump_stats('cluster_list.pstat')
    return result


def test_r_tree(tau, dim=2, sample_size=2000,sq_metric=mtr.sq_dist_2, profiler=False, build_method='str', leaf_points=40, inner_children=20):
    print("---- Testing R-Tree range search for %s sample points" % sample_size, " and tau = ", tau, " --- ")
    start_time = time.time()
    r_tree:RTree = None
    if build_method == 'str':
        r_tree = bulk_load_r_tree_str(point_list, leaf_points, leaf_points//2, inner_children, inner_children//2, sq_metric)
    if build_method == 'morton':
        r_tree = bulk_load_r_tree_morton(point_list, leaf_points, leaf_points//2, inner_children, inner_children//2, sq_metric)
    if build_method == 'single':
        r_tree:RTree = RTree(None, inner_children, inner_children//2, leaf_points, leaf_points//2, sq_metric)
        for p in point_list:
            r_tree.insert(p)
    print("--- %s seconds " % (time.time() - start_time), " for building R-Tree ---")
    print("Depth", r_tree.get_depth())
    print("Correct", r_tree.is_correct())
    print("Volume:", r_tree.get_volume())
    print("Average Leaf volume, Leaf number:", r_tree.average_leaf_volume())
    print("Number of leafs")
    result = {}
    smpl = random.sample(point_list, sample_size)
    start_time = time.time()
    pr = profile.Profile()
    if profiler:
        pr.enable()
    for p in smpl:
        result[p] = r_tree.get_points_in_rect_it(Rect((p[0]-tau, p[1]-tau), (p[0]+tau, p[1]+tau)))
    if profiler:
        pr.disable()
    print("--- %s seconds " % ((time.time() - start_time)/sample_size), " per point query")
    if profiler:
        pr.dump_stats('r_tree.pstat')
    return result


def test_sequential(tau, dim=2, sample_size=2000,sq_metric=mtr.sq_dist_2, profiler=False):
    print("---- Testing sequential range search for %s sample points" % sample_size, " and tau = ", tau, " --- ")
    start_time = time.time()
    result = {}
    smpl = random.sample(point_list, sample_size)
    start_time = time.time()
    pr = profile.Profile()
    if profiler:
        pr.enable()
    tau_sq = tau*tau
    for p in smpl:
        result[p] = [y for y in point_list if sq_metric(y,p) <= tau_sq]
    if profiler:
        pr.disable()
    print("--- %s seconds " % ((time.time() - start_time)/sample_size), " per point query")
    if profiler:
        pr.dump_stats('seq.pstat')
    return result

import scipy.spatial.kdtree
def test_scipy_kd_tree(tau, dim=2, sample_size=2000, profiler=False, leaf_points=100):
    print("---- Testing Scipy KD-Tree range search for %s sample points" % sample_size, " and tau = ", tau, " --- ")
    start_time = time.time()
    kd_tree = scipy.spatial.kdtree.KDTree(point_list, leafsize=leaf_points)
    print("--- %s seconds " % (time.time() - start_time), " for building KD-Tree ---")
    result = {}
    smpl = random.sample(point_list, sample_size)
    start_time = time.time()
    pr = profile.Profile()
    if profiler:
        pr.enable()
    result = list(map(lambda x:list(map(lambda y:point_list[y],x)),kd_tree.query_ball_tree(kd_tree, tau)))
    #for p in smpl:
    #    result[p] = list(map(lambda x:point_list[x], kd_tree.query_ball_point(p, tau)))
    if profiler:
        pr.disable()
    print("--- %s seconds " % ((time.time() - start_time)/sample_size), " per point query")
    if profiler:
        pr.dump_stats('kd_scipy_tree.pstat')
    return result


def test_correctness_of_range_query(result, tau, sample_size=10, metric=mtr.dist_mem, sq_metric=mtr.sq_dist):
    smpl = random.sample(result.keys(), sample_size)
    tau_sq = tau*tau
    for p in smpl:
        l = [y for y in point_list if sq_metric(y,p) <= tau_sq]
        if set(l) != set(result[p]):
            print("Not correct at ", p)
            print("True list is subset of computed list: ", set(l).issubset(set(result[p])))
            print("Computed list is subset of true list: ", set(result[p]).issubset(set(l)))
            print("Computed list length | True list length: ", result[p].__len__(), " | ", l.__len__())

            return
    print("Correct")


def compare_results(resutl_1, result_2):
    if result_1.__len__() != result_2.__len__():
        print("Results are not of equal length:", result_1.__len__(), result_2.__len__())
        return
    for p in result_1:
        if set(result_1[p]) != set(result_2[p]):
            print("Not correct at ", p)
            print("True list is subset of computed list: ", set(result_1[p]).issubset(set(result_2[p])))
            print("Computed list is subset of true list: ", set(result_2[p]).issubset(set(resutl_1[p])))
            print("result_1 length | result_2 length: ", result_1[p].__len__(), " | ", result_2[p].__len__())

            return
    print("Equal")

tau = 0.05

result = test_scipy_kd_tree(tau, sample_size=point_list.__len__(), profiler=False, leaf_points=100)
#result_1 = test_m_tree(tau,sample_size=100, dim=dimension, profiler=False,sq_metric=use_sq_metric, leaf_points=20, inner_children=7)
#result = test_sequential(tau, sample_size=1000)
#result = test_cluster_list(tau, profiler=False,disk_radius=0.3)
#result_2 = test_r_tree(tau,sample_size=100, dim=dimension, build_method='str', sq_metric=use_sq_metric, leaf_points=20, inner_children=10, profiler=True)
test_correctness_of_range_query(result, tau, sample_size=10)
#compare_results(result_1, result_2)
'''
m_tree = MTree(mtr.dist_mem, 6, 20, sq_metric=mtr.sq_dist_2)
for p in point_list:
    m_tree.insert(p)

plt.axis([0, 1, 0, 1])
ax = plt.gca()
ax.set_autoscale_on(False)
queue = deque([m_tree.root])
while queue:
    node = queue.pop()
    if node.is_leaf:
        ax.add_artist(plt.Circle(node.obj, node.covering_radius, color='r', fill=False))
    else:
        queue.extend(node.subtree)

plt.scatter([x[0] for x in point_list], [x[1] for x in point_list], color='b')
plt.savefig("m_tree_no_morton_bananas-1-2d.png")
'''