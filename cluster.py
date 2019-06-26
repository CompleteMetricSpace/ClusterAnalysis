import csv, os
import math
from collections import deque
import time
from m_tree import MTree
from r_tree import RTree, bulk_load_r_tree_str
from kd_tree import KDTree, create_kd_tree
import metrics as mtr
import disjoint_set as dst
from morton import morton_code
import logging
import multiprocessing
from multiprocessing import Pool


dataListGlobal = []
index_inverse_global = {}

def neighbour(t):
    tree, rg, tau, type = t
    nd = []
    if type == 'm':
        for p in rg:
            nd.append(tree.get_points_in_range_it(p, tau))
    else:
        for p in rg:
            nd.append([index_inverse_global[q] for q in tree.get_points_in_range_it(dataListGlobal[p], tau)])
    return nd


def chunk(num, k):
    av = int(math.ceil(num / k))
    l = []
    for i in range(0,num, av):
        l.append(range(i,min(i+av,num)))
    return l


def get_filename(name):
    head, tail = os.path.split(name)
    if not tail.endswith(".train.csv"):
        raise Exception("Filename must be of the form 'name.train.csv' if output name is not specified.")
    return tail[0:tail.__len__()-10]


def  compute_neighbours(number_of_threads, n, tree, tree_structure, tau):
    if number_of_threads > 1:
        pl = Pool()
        index_list = chunk(n, number_of_threads)
        input_list = [(tree, ind, tau, tree_structure) for ind in index_list]
        output = pl.map(neighbour, input_list)
        # neighbour_dict = m_tree.get_points_in_range_it_list(dataList, tau)
        neighbour_radius_dict = []
        for l in output:
            neighbour_radius_dict.extend(l)
    else:
        neighbour_radius_dict = neighbour((tree, range(n), tau, tree_structure))
    return neighbour_radius_dict



# The main algorithm

def cluster(name, epsilon: float, delta: float, tau: float, tree_structure='automatic', threads=1,
            sequential_search_threshold=2000, return_list=False, epsilon_multiple=None, output_name=None):
    """
    The main clustering algorithm "Adaptive Density Level Set Clustering"
    :param name: the file name of the csv file
    :param epsilon: a number > 0
    :param delta: a number > 0
    :param tau: a number > 0
    :param tree_structure: a tree structure type 'm', 'r', 'kd' or 'automatic'
    :param threads: the number of threads > 0 or 'automatic'
    :param sequential_search_threshold: the threshold below which sequential search for computing density will be applied
    :param return_list: true or false
    :return: if return_list is true, then this returns a list of components, None otherwise
    """
    global dataListGlobal, index_inverse_global

    if output_name is None:
        output_name = get_filename(name)

    neighbour_radius = max(delta, tau)

    # Choose tree structure if tree structure is automatic
    if tree_structure == 'automatic':
        if neighbour_radius > 0.05:
            tree_structure = 'm'
        else:
            if neighbour_radius > 0.01:
                tree_structure = 'kd'
            else:
                tree_structure = 'r'

    if tree_structure not in ['r', 'm', 'kd']:
        raise Exception("Tree structure can only be m, r, kd")

    logging.info("Tree structure chosen: %s", tree_structure)

    # Setup
    dataList = []
    with open(name) as csvfile:
        data = csv.reader(csvfile)
        for row in data:
            dataList.append(tuple([float(s) for s in row]))


    # Only a portion of the data is taken
    portion = dataList.__len__()
    #portion = 500
    dataList = dataList[0:portion]

    original_data = dataList.copy()


    if tree_structure == 'm':
        stime = time.time()
        logging.info("Sorting using z-curve...")
        # Sort dataList:
        dataList.sort(key=lambda x: morton_code(x))
        logging.info("Data set sorted. Time for sorting: %f", (time.time() - stime))

    #Define constants
    n: int = len(dataList)
    d: int = len(dataList[0])
    dimension_const: float = n * delta ** d

    dataListGlobal = dataList

    mtr.dataList = dataList

    #Define sq_metric
    sq_metric_index = mtr.choose_sq_metric_index(d)
    sq_metric = mtr.choose_sq_metric(d)

    #Define metric
    metric_index = mtr.choose_metric_index(d)



    # Create index inverse
    index_inverse_global = {dataList[i]: i for i in range(n)}


    # Build tree
    stime = time.time()
    if tree_structure == 'm':
        logging.info("Filling M-Tree...")
        tree = MTree(metric_index, 6, 20, sq_metric=sq_metric_index)
        for p in range(n):
            tree.insert(p)
        logging.info("Done filling M-Tree. Time for filling: %f",  (time.time()-stime))
    if tree_structure == 'r':
        logging.info("Filling R-Tree...")
        #Make a copy, since bulk_load_r_tree sorts dataList
        dataListCopy = dataList.copy()
        tree = bulk_load_r_tree_str(dataListCopy, 40, 20, 20, 10, sq_metric=sq_metric)
        logging.info("Done filling R-Tree. Time for filling: %f", (time.time() - stime))
    if tree_structure == 'kd':
        logging.info("Constructing KD-Tree...")
        tree = create_kd_tree(dataList, 50, sq_metric=sq_metric)
        logging.info("KD-Tree constructed. Time for construction: %f", (time.time() - stime))


    stime = time.time()
    logging.info("Precomputing neighbourhood list...")
    number_of_cores = multiprocessing.cpu_count()
    if threads == 'automatic':
        if n > 5000:
            number_of_threads = math.ceil(number_of_cores/2)
        else:
            if n > 2500:
                number_of_threads = min(number_of_cores, 2)
            else:
                number_of_threads = 1
    else:
        number_of_threads = threads

    logging.info("Number of cores detected: %s. Using %s thread(s)", number_of_cores, number_of_threads)


    neighbour_radius_dict = compute_neighbours(number_of_threads, n, tree, tree_structure, tau)


    logging.info("max(tau,delta)-neighbourhood list computed. Time for computation: %f",  (time.time()-stime))
    average_neighbour_number  = sum([neighbour_radius_dict[p].__len__() for p in range(n)])/n
    logging.info("Average number of neighbours: %d", average_neighbour_number)
    stime = time.time()
    if math.isclose(delta,neighbour_radius) and math.isclose(tau, neighbour_radius):
        logging.info("Delta = Tau, saving distance calculations")
        neighbour_dict = neighbour_radius_dict
        density_dict = [neighbour_radius_dict[i].__len__() for i in range(n)]
    else:
        if math.isclose(delta, neighbour_radius):
            logging.info("Tau is smaller than delta, computing tau-neighbours...")
            if average_neighbour_number <= sequential_search_threshold:
                #Compute tau-neighbours sequentially
                neighbour_dict = [[] for x in range(n)]
                tau_sq = tau*tau
                for p in range(n):
                    for y in neighbour_radius_dict[p]:
                        if y <= p and sq_metric_index(y, p) <= tau_sq:
                            neighbour_dict[p].append(y)
                            neighbour_dict[y].append(p)
            else:
                #Compute tau-neighbours through the tree structure
                neighbour_dict = compute_neighbours(number_of_threads, n, tree, tree_structure, tau)
            density_dict = [l.__len__() for l in neighbour_radius_dict]
            logging.info("Tau-neighbours computed. Time for computation: %f", (time.time()-stime))
        else:
            logging.info("Delta is smaller than tau, computing density...")
            if average_neighbour_number <= sequential_search_threshold:
                #Compute delta-neighbours sequentially
                density_dict = [0 for x in range(n)]
                delta_sq = delta * delta
                for p in range(n):
                    for y in neighbour_radius_dict[p]:
                        if y <= p and sq_metric_index(y, p) <= delta_sq:
                            density_dict[p] += 1
                            density_dict[y] += 1
            else:
                # Compute delta-neighbours through the tree structure
                d_dict = compute_neighbours(number_of_threads, n, tree, tree_structure, delta)
                density_dict = [d_dict[i].__len__() for i in range(n)]
            neighbour_dict = neighbour_radius_dict
            logging.info("Density computed. Time for computation: %f",  (time.time()-stime))


    #Allocate a disjoint-set (find-union) data structure
    disjoint_set_dict = [None]*n
    rank_dict = [None]*n

    for x in range(n):
        #Initialize all disjoint sets
        disjoint_set_dict[x] = x
        rank_dict[x] = 1


    inc_constant: float = 1 / dimension_const

    upper_bound_on_density = max([density_dict[p] for p in range(n)])/dimension_const
    logging.info("Upper bound on density: %f. Required iterations: %f", upper_bound_on_density, upper_bound_on_density*dimension_const)

    if epsilon_multiple != None:
        epsilon = math.sqrt(upper_bound_on_density/dimension_const)*epsilon_multiple

    rho: float = upper_bound_on_density

    old_level_set = {}

    number_of_connected_components = 0

    #Save the last branching point

    saved_rho = 0
    saved_disjoint_set_dict = []
    saved_level_set_rho = {}
    saved_number_of_components = 0
    logging.info("Entering loop")
    start_time = time.time()
    # The incrementation is decreased until rho = 0
    while rho >= -1.0e-9:
        logging.info("--- new loop ---")
        epsrho: float = rho + 2 * epsilon
        logging.info("Computing level set...")

        level_set_rho = {y:True for y in range(n) if density_dict[y]/dimension_const >= rho}
        if rho == upper_bound_on_density:
            number_of_connected_components = level_set_rho.__len__()

        logging.info("Level set computed")

        logging.info("Computing connected components...")

        #For each new element call union
        for y in level_set_rho:
            if y not in old_level_set:
                number_of_connected_components += 1
                for z in neighbour_dict[y]:
                    if z in level_set_rho:
                        if dst.union(y, z, disjoint_set_dict, rank_dict):
                            number_of_connected_components -= 1

        logging.info("Connected components computed: %d", number_of_connected_components)

        logging.info("Removing insignificant components...")

        #TODO: Optimize this part, i.e. persisting_components should be an integer
        persisting_components = {}
        for y in level_set_rho:
            if density_dict[y]/dimension_const >= epsrho:
                persisting_components[dst.get_root(y, disjoint_set_dict)] = True

        logging.info("Insignificant components removed. Now there are: %d", persisting_components.__len__())

        if persisting_components.__len__() > 1:
            logging.info("Branching point, data saved.")
            #Save the branching point
            saved_rho = rho
            saved_disjoint_set_dict = disjoint_set_dict.copy()
            saved_level_set_rho = level_set_rho.copy()
            saved_number_of_components = persisting_components.__len__()

        # Decrease rho
        rho -= inc_constant
        old_level_set = level_set_rho
        logging.info("Decreasing rho to: %f", rho)

    logging.info("--- %f seconds in loop ---", (time.time()-start_time))

    logging.info("Components detected %s", saved_number_of_components)
    logging.info("Writing to file...")

    #Determine component dictionary

    component_dict = [None]*n
    counter = 1
    for y in range(n):
        if y in saved_level_set_rho:
            root = dst.get_root(y, saved_disjoint_set_dict)
            if not component_dict[root] is None:
                component_dict[y] = component_dict[root]
            else:
                component_dict[root] = counter
                component_dict[y] = counter
                counter += 1
        else:
                component_dict[y] = 0

    #Write to file in the original order
    failed = saved_level_set_rho.__len__() == 0

    with open(output_name+'.result.csv', 'w') as file:
        for p in original_data:
            if failed:
                file.write(str(1))
            else:
                file.write(str(component_dict[index_inverse_global[p]]))
            file.write(', ')
            file.write(', '.join(map(str,p)))
            file.write('\n')
    file.close()

    #Create list of components
    if return_list:
        list_of_components = [[] for i in range(counter)]
        for y in range(n):
            list_of_components[component_dict[y]].append(dataList[y])
        return list_of_components

    logging.info("Done.")


logging.getLogger().setLevel(logging.INFO)
start_time = time.time()

list_of_components = cluster(
    name = '/home/andrei/Documents/University/Computerpraktikum/Python/ExtractedFiles/Cluster/Artificial/bananas-5-2d.csv',
    delta=0.045, epsilon=0.49, tau=0.09, threads=2, output_name='test', return_list=True)
#For banana-5-2d take delta = 0.045, tau = 0.09, epsilon = 0.49
print("--- %f seconds total ---" % (time.time() - start_time))
import plots
plots.plot_components(list_of_components)
