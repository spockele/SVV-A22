"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as axes3d
import numpy as np


def read_inp():
    """
    Reads the Abaqus input file
    :return:
    """
    path = "B737.inp"
    with open(path) as f:
        nodes = False
        data_opt = {}

        for line in f:
            if line == "*Element, type=S4R\n":
                break

            if nodes:
                data_line = line.strip("\n").replace(' ', '')
                data_line = data_line.split(",")
                for i, data_pt in enumerate(data_line):
                    data_line[i] = float(data_pt) if i > 0 else int(data_pt)

                data_opt[data_line[0]] = data_line[1:]

            if line == "*Node\n":
                nodes = True

    return data_opt


def plot_aileron(node_dict: dict):
    """
    Plot the node coordinates of the inp file
    :param node_dict: the return from read_inp()
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xs = np.zeros((len(node_dict),))
    ys = np.zeros((len(node_dict),))
    zs = np.zeros((len(node_dict),))

    for node_num, node_coord in node_dict.items():
        xs[node_num - 1] = node_coord[0]
        ys[node_num - 1] = node_coord[1]
        zs[node_num - 1] = node_coord[2]

    ax.scatter(xs, ys, zs)

    plt.show()


if __name__ == '__main__':
    node = read_inp()
    plot_aileron(node)
