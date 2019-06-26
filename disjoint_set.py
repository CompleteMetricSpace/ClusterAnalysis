
#The disjoint set data structure
#See wikipedia for more information

def get_root(x, disjoint_set_dict):
    y = x
    while disjoint_set_dict[y] != y:
        #y, disjoint_set_dict[y] = disjoint_set_dict[y], disjoint_set_dict[disjoint_set_dict[y]]
        y = disjoint_set_dict[y]
    #disjoint_set_dict[x] = y
    return y


def union(x, y, disjoint_set_dict, rank_dict):
    root_x = get_root(x, disjoint_set_dict)
    root_y = get_root(y, disjoint_set_dict)
    if root_x == root_y:
        return False
    if rank_dict[root_x] > rank_dict[root_y]:
        disjoint_set_dict[root_y] = root_x
    else:
        disjoint_set_dict[root_x] = root_y
        if (rank_dict[root_x] == rank_dict[root_y]):
            rank_dict[root_y] = rank_dict[root_y] + 1
    return True