"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as axes3d
import numpy as np


class Element():
    def __init__(self, num, nodes, node_dct):
        self.num = num
        self.nodes = nodes

        self.corners = []
        for node in nodes:
            self.corners.append(node_dct[node])

        self.corners = self.corners[:-2] + [self.corners[-1]] + [self.corners[-2]]

        self.xs, self.ys, self.zs = np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))
        for i, corner in enumerate(self.corners):
            self.xs[i % 2, i // 2], self.ys[i % 2, i // 2], self.zs[i % 2, i // 2] = corner[0], corner[1], corner[2]

    def __repr__(self):
        return f'Element: nodes: {self.nodes}\ncorners: {self.corners}\n'

    def __str__(self):
        return f'Element: nodes: {self.nodes}, corners: {self.corners}'

    def plot_self(self, ax):
        ax.plot_wireframe(self.xs, self.zs, self.ys)


def read_inp():
    """
    Reads the Abaqus input file
    :return:
    """
    path = "B737.inp"
    with open(path) as f:
        nodes = False
        elements = False
        s = False
        node_dct = {}
        elem_lst = []

        for line in f:
            if line == "*Nset, nset=Skin\n":
                elements = False
                s = False
                break

            if line == "*Element, type=S4R\n":
                nodes = False
                elements = True
                s = False

            if line == "*Node\n":
                nodes = True
                s = False

            if elements and s:
                data_line = line.strip("\n").replace(' ', '')
                data_line = data_line.split(',')
                for i, data_pt in enumerate(data_line):
                    data_line[i] = int(data_pt)

                elem_lst.append(Element(data_line[0], data_line[1:], node_dct))

            if nodes and s:
                data_line = line.strip("\n").replace(' ', '')
                data_line = data_line.split(",")
                for i, data_pt in enumerate(data_line):
                    data_line[i] = float(data_pt) if i > 0 else int(data_pt)

                node_dct[data_line[0]] = tuple(data_line[1:])

            s = True

    return node_dct, elem_lst


def read_rpt():
    path = "B737_working_copy.rpt"
    with open(path) as f:
        case, tpe = False, False
        for line in f:
            if line == "\n":
                case, tpe = False, False

            if case and tpe:
                print(line)

            if "Step" in line:
                case, tpe = line.strip("Step: ").split(", ")

    return None, None, None


if __name__ == '__main__':
    case1, case2, case3 = read_rpt()
    #node, elem = read_inp()
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111, projection='3d')
    #for el in elem:
    #    el.plot_self(ax1)

    #plt.show()
