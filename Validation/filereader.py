"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import matplotlib.colors as color
import matplotlib.cm as cm
import mpl_toolkits.mplot3d as axes3d
import numpy as np


class Element:
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


class Node:
    def __init__(self, num, xyz, data_dict, rf, s=None):
        self.num = num
        self.x, self.y, self.z = xyz
        if s is None:
            s = "displacement"

        self.U_Bend = data_dict["Bending"][s][self.num]
        self.U_JBen = data_dict["Jam_Bent"][s][self.num]
        self.U_JStr = data_dict["Jam_Straight"][s][self.num]

        self.xd, self.yd, self.zd = None, None, None

        if rf:
            self.R_Bend = data_dict["Bending"]["reaction"][self.num]
            self.R_JBen = data_dict["Jam_Bent"]["reaction"][self.num]
            self.R_JStr = data_dict["Jam_Straight"]["reaction"][self.num]

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

    def plot_self(self, ax, case, u=False, color="red"):
        if u:
            self.displace(case)
            ax.scatter([self.xd], [self.zd], [self.yd], color=color)
        else:
            ax.scatter([self.x], [self.z], [self.y], color=color)

    def plot_quiver(self, ax, case, u=False, sc=1):
        if case == "Bending":
            _, rx, ry, rz = self.R_Bend
        elif case == "Jam_Bent":
            _, rx, ry, rz = self.R_JBen
        elif case == "Jam_Straight":
            _, rx, ry, rz = self.R_JStr

        if u:
            self.displace(case)
            ax.scatter([self.xd], [self.zd], [self.yd], color='red')
            ax.quiver([self.xd], [self.zd], [self.yd], rx*sc, rz*sc, ry*sc, color="red")
        else:
            ax.scatter([self.x], [self.z], [self.y], color='red')
            ax.quiver([self.x], [self.z], [self.y], rx*sc, rz*sc, ry*sc, color="red")


def read_inp(data_dict):
    """
    Reads the Abaqus input file
    :return:
    """
    path = "B737_working_copy.inp"
    with open(path) as f:
        nodes = False
        elements = False
        assembly = False
        s = False
        node_dct = {}
        elem_lst = []
        asem_lst = []

        for line in f:
            if line == "** ASSEMBLY\n":
                elements = False
                assembly = True
                s = False

            if line == "*Element, type=S4R\n":
                nodes = False
                elements = True
                s = False

            if line == "*Node\n":
                nodes = True
                s = False

            if assembly and s:
                data_line = line.strip("\n").replace(' ', '')
                data_line = data_line.split(',')
                for i, data_pt in enumerate(data_line):
                    data_line[i] = float(data_pt) if i > 0 else int(data_pt)

                asem_lst.append(Node(data_line[0], tuple(data_line[1:]), data_dict, True, "ASSEMBLY"))

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

                node_dct[data_line[0]] = Node(data_line[0], tuple(data_line[1:]), data_dict, False)

            s = True

    return node_dct, elem_lst, asem_lst


def read_rpt():
    path = "B737_working_copy.rpt"
    with open(path) as f:
        case, tpe = False, False
        data_dict = {
                     "Bending":      {"VM_S12": {}, "displacement": {}, "reaction": {}, "ASSEMBLY": {}},
                     "Jam_Bent":     {"VM_S12": {}, "displacement": {}, "reaction": {}, "ASSEMBLY": {}},
                     "Jam_Straight": {"VM_S12": {}, "displacement": {}, "reaction": {}, "ASSEMBLY": {}}
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
    node, elem, asem = read_inp(data)

    mm = VMi[case] if tpe == "VMi" else S12[case]

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.set_ylim(0-600, 100)
    ax1.set_zlim(0-350, 350)

    #for el in elem:
    #    el.plot_self(ax1, case, tpe, VMi, S12)

    for no in asem:
        no.plot_self(ax1, case, u=True)
        scale = 100 if no.num in range(5, 16) else 1
        no.plot_quiver(ax1, case, u=True, sc=scale)

    for n, no in node.items():
        if n in (11,   14,   39,   40,   41,   42,  207,  208,  209,  210,  211,  918,  919,  920,  921,  922,
                 923,  924,  925,  926,  927,  928,  929,  930,  931,  932,  933,  934,  935,  955,  963,  964,
                 965,  966,  967,  968,  969,  970,  971,  972,  973,  974,  975,  976,  977,  978,  979,  980,
                 981,  982,  983,  984,  985,  986,  987, 1028, 1029, 1030, 1031, 1032, 1033,
                 1, 4, 18, 19, 24, 26, 52, 53, 54, 55, 56, 305, 306, 307, 308, 309,
                 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 458, 461, 462,
                 463, 464, 465, 532, 533, 534, 535, 536, 537, 538, 541, 566, 567, 568, 569, 570,
                 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583,
                 7, 8, 27, 29, 30, 31, 162, 163, 164, 165, 166, 553, 554, 555, 556, 557,
                 565, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622,
                 623, 624, 625, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685,
                 686, 687, 688, 689, 690, 721, 723, 724, 725, 726, 727, 728, 729
                 ):
            no.plot_self(ax1, case, u=True, color="green")

        if n in (12,  13,  15,  36,  37,  38, 196, 197, 198, 199, 200, 246, 247, 248, 249, 250,
                 823, 876, 877, 878, 879, 880, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890,
                 891, 892, 893, 894, 895, 896, 897, 898, 899, 900, 901, 902, 903, 904, 905, 906,
                 907, 908, 909, 910, 911, 948, 956, 957, 958, 959, 960, 961, 962,
                 9, 10, 32, 33, 34, 35, 179, 180, 181, 182, 183, 697, 698, 699, 700, 701,
                 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 722, 730, 731,
                 732, 733, 734, 735, 736, 749, 750, 751, 752, 753, 788, 824, 825, 826, 827, 828,
                 829, 830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841,
                 5, 6, 21, 22, 25, 28, 104, 105, 106, 107, 108, 390, 391, 392, 393, 394,
                 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 460, 478, 479,
                 480, 481, 482, 483, 484, 542, 543, 544, 545, 546, 558, 584, 585, 586, 587, 588,
                 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601,
                 2, 3, 16, 17, 20, 23, 45, 46, 47, 48, 49, 285, 286, 287, 288, 289,
                 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 325, 326, 327,
                 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 455,
                 459, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477
                 ):
            no.plot_self(ax1, case, u=True, color="blue")

    n = color.Normalize(vmin=mm[0], vmax=mm[1])
    m = cm.ScalarMappable(norm=n, cmap=plt.viridis())
    plt.colorbar(m)

    ax1.set_xlabel("X")
    ax1.set_ylabel("Z")
    ax1.set_zlabel("Y")
    plt.show()
