import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import colorsys

def get_distinct_colors(k):
    hsvTuples = [(x * 1.0 / k, 0.5, 0.5) for x in range(k)]
    return list(map(lambda x: colorsys.hsv_to_rgb(*x), hsvTuples))


def plot_components(list_of_components, coord=(0,1), file=None):
    """
    Plots the components
    :param list_of_components: a list of lists of points
    :param coord: a coordinate vector, i.e. if coord = (a,b) then points of the form (x[a],x[b]) are plotted
    :return: None
    """
    k = list_of_components.__len__()
    colorList = get_distinct_colors(k)
    plt.axis([0, 1, 0, 1])
    ax = plt.gca()
    ax.set_autoscale_on(False)
    counter = 0
    for component in list_of_components:
        plt.scatter([x[coord[0]] for x in component], [x[coord[1]] for x in component], color=colorList[counter])
        counter +=1
    color_patch = mpatches.Patch(color=colorList[0], label='Non classified points')
    plt.legend(handles=[color_patch])
    if file is None:
        plt.show()
    else:
        plt.savefig(file)

