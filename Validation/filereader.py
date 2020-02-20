"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import matplotlib.colors as color
import matplotlib.cm as cm
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
        self.xs, self.ys, self.zs = None, None, None

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

    def plot_self(self, ax, case, tpe, VMi_mm, S12_mm):
        self.xs, self.ys, self.zs = np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))
        for i, corner in enumerate(self.corners):
            corner.displace(case)
            self.xs[i % 2, i // 2], self.ys[i % 2, i // 2], self.zs[i % 2, i // 2] = corner.xd, corner.yd, corner.zd

        if case == "Bending":
            col = (self.VMi_Bend - VMi_mm[case][0]) / (VMi_mm[case][1] - VMi_mm[case][0]) if tpe == "VMi" else (self.S12_Bend - S12_mm[case][0]) / (S12_mm[case][1] - S12_mm[case][0])
        elif case == "Jam_Bent":
            col = (self.VMi_JBen - VMi_mm[case][0]) / (VMi_mm[case][1] - VMi_mm[case][0]) if tpe == "VMi" else (self.S12_JBen - S12_mm[case][0]) / (S12_mm[case][1] - S12_mm[case][0])
        elif case == "Jam_Straight":
            col = (self.VMi_JStr - VMi_mm[case][0]) / (VMi_mm[case][1] - VMi_mm[case][0]) if tpe == "VMi" else (self.S12_JStr - S12_mm[case][0]) / (S12_mm[case][1] - S12_mm[case][0])

        #print(col)
        cmap = plt.get_cmap('viridis')
        col = [[cmap(col)]]

        ax.plot_surface(self.xs, self.zs, self.ys, facecolors=col)


class Node():
    def __init__(self, num, xyz, data_dict):
        self.num = num
        self.x, self.y, self.z = xyz 
        self.U_Bend = data_dict["Bending"]["displacement"][self.num]
        self.U_JBen = data_dict["Jam_Bent"]["displacement"][self.num]
        self.U_JStr = data_dict["Jam_Straight"]["displacement"][self.num]

        self.xd, self.yd, self.zd = None, None, None

    def __repr__(self):
        return f'Node {self.num}: ({self.x}, {self.y}, {self.z}), U: {self.U_Bend}, {self.U_JBen}, {self.U_JStr}\n'

    def displace(self, case):
        if case == "Bending":
            _, dx, dy, dz = self.U_Bend
        elif case == "Jam_Bent":
            _, dx, dy, dz = self.U_JBen
        elif case == "Jam_Straight":
            _, dx, dy, dz = self.U_JStr
        else:
            dx, dy, dz = 0, 0, 0

        self.xd, self.yd, self.zd = self.x + dx, self.y + dy, self.z + dz


    def plot_self(self, ax, u=False):
        if u:
            ax.scatter([self.xd], [self.zd], [self.yd], color='red')
        else:
            ax.scatter([self.x], [self.z], [self.y], color='red')


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

                node_dct[data_line[0]] = Node(data_line[0], tuple(data_line[1:]), data_dict)

            s = True

    return node_dct, elem_lst


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

            if "Step" in line:
                case, tpe = line.strip("\n").replace("Step: ", "").split(", ")

        VMi_lst = {}
        S12_lst = {}
        for case in data_dict:
            VMi_lst[case] = [(data_dict[case]["VM_S12"][node][0] + data_dict[case]["VM_S12"][node][1]) / 2 for node in data_dict[case]["VM_S12"]]
            S12_lst[case] = [(data_dict[case]["VM_S12"][node][2] + data_dict[case]["VM_S12"][node][3]) / 2 for node in data_dict[case]["VM_S12"]]

        VMi = {case: (min(VMi_lst[case]), max(VMi_lst[case])) for case in VMi_lst}
        S12 = {case: (min(S12_lst[case]), max(S12_lst[case])) for case in S12_lst}

    return data_dict, VMi, S12


if __name__ == '__main__':
    ipt = input("Load case and stress to display: ")
    case, tpe = ipt.split(", ")
    data, VMi, S12 = read_rpt()
    node, elem = read_inp(data)

    mm = VMi[case] if tpe == "VMi" else S12[case]

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_ylim(0-600, 100)
    ax1.set_zlim(0-350, 350)

    for el in elem:
        el.plot_self(ax1, case, tpe, VMi, S12)

    #for no in node:
    #    no.plot_self(ax1)

    n = color.Normalize(vmin=mm[0], vmax=mm[1])
    m = cm.ScalarMappable(norm=n, cmap=plt.viridis())
    plt.colorbar(m)
    plt.show()
