"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as axes3d
import numpy as np


class Element():
    def __init__(self, num, nodes, node_dct, data_dict):
        self.num = num
        self.nodes = nodes

        self.corners = []
        for node in nodes:
            self.corners.append(node_dct[node])

        self.corners = self.corners[:-2] + [self.corners[-1]] + [self.corners[-2]]

        self.xs, self.ys, self.zs = np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))
        for i, corner in enumerate(self.corners):
            self.xs[i % 2, i // 2], self.ys[i % 2, i // 2], self.zs[i % 2, i // 2] = corner[0], corner[1], corner[2]

        self.VMi_Bend = (data_dict["Bending"]["VM_S12"][self.num][0] + data_dict["Bending"]["VM_S12"][self.num][1]) / 2
        self.S12_Bend = (data_dict["Bending"]["VM_S12"][self.num][2] + data_dict["Bending"]["VM_S12"][self.num][3]) / 2

        self.VMi_JBen = (data_dict["Jam_Bent"]["VM_S12"][self.num][0] + data_dict["Jam_Bent"]["VM_S12"][self.num][1]) / 2
        self.S12_JBen = (data_dict["Jam_Bent"]["VM_S12"][self.num][2] + data_dict["Jam_Bent"]["VM_S12"][self.num][3]) / 2

        self.VMi_JStr = (data_dict["Jam_Straight"]["VM_S12"][self.num][0] + data_dict["Jam_Straight"]["VM_S12"][self.num][1]) / 2
        self.S12_JStr = (data_dict["Jam_Straight"]["VM_S12"][self.num][2] + data_dict["Jam_Straight"]["VM_S12"][self.num][3]) / 2

    def __repr__(self):
        return f'Element {self.num}: nodes: {self.nodes},'\
               f' VMi: ({self.VMi_Bend}, {self.VMi_JBen}, {self.VMi_JStr}),'\
               f' S12: ({self.S12_Bend}, {self.S12_JBen}, {self.S12_JStr})\n'

    def plot_self(self, ax):
        ax.plot_wireframe(self.xs, self.zs, self.ys)


class Node():
    def __init__(self, num, xyz, data_dict):
        self.num = num
        self.x, self.y, self.z = xyz
        self.U_Bend = data_dict["Bending"]["displacement"][self.num]
        self.U_JBen = data_dict["Jam_Bent"]["displacement"][self.num]
        self.U_JStr = data_dict["Jam_Straight"]["displacement"][self.num]

    def __repr__(self):
        return f'Node {self.num}: ({self.x}, {self.y}, {self.z}), U: {self.U_Bend}, {self.U_JBen}, {self.U_JStr}\n'


def read_inp(data_dict):
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
        node_lst = []
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

                elem_lst.append(Element(data_line[0], data_line[1:], node_dct, data_dict))

            if nodes and s:
                data_line = line.strip("\n").replace(' ', '')
                data_line = data_line.split(",")
                for i, data_pt in enumerate(data_line):
                    data_line[i] = float(data_pt) if i > 0 else int(data_pt)

                node_lst.append(Node(data_line[0], tuple(data_line[1:]), data_dict))
                node_dct[data_line[0]] = tuple(data_line[1:])

            s = True

    return node_lst, elem_lst


def read_rpt():
    path = "B737_working_copy.rpt"
    with open(path) as f:
        case, tpe = False, False
        data_dict = {
                     "Bending":      {"VM_S12": {}, "displacement": {}, "reaction": {}},
                     "Jam_Bent":     {"VM_S12": {}, "displacement": {}, "reaction": {}},
                     "Jam_Straight": {"VM_S12": {}, "displacement": {}, "reaction": {}}
                     }
        for line in f:
            if line == "\n":
                case, tpe = False, False

            if case and tpe:
                data_line = line.strip("\n").split(",")
                data_line = [int(data_line[0])] + [float(data) for data in data_line[1:]]
                data_dict[case][tpe][data_line[0]] = tuple(data_line[1:])
                print(case, ':', tpe, ':', data_line[0], ':', data_dict[case][tpe][data_line[0]])

            if "Step" in line:
                case, tpe = line.strip("\n").replace("Step: ", "").split(", ")

    return data_dict


if __name__ == '__main__':
    data = read_rpt()
    node, elem = read_inp(data)
    print(node)
    print(elem)
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_xlim(0, 1200)
    ax1.set_ylim(0-600, 600)
    ax1.set_zlim(0-600, 600)
    for el in elem:
        el.plot_self(ax1)

    plt.show()
