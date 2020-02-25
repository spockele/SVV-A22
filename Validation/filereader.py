"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import matplotlib.colors as color
import matplotlib.cm as cm
import mpl_toolkits.mplot3d as axes3d
import numpy as np


class Element:
    def __init__(self, num, nodes, node_dct, data_dict, VMi_mm, S12_mm):
        self.num = num
        self.nodes = nodes

        self.corners = []
        for n in nodes:
            self.corners.append(node_dct[n])

        self.corners = self.corners[:-2] + [self.corners[-1]] + [self.corners[-2]]
        self.xs, self.ys, self.zs = None, None, None
        self.z, self.y, self.x = None, None, None

        self.data = {case: {"VMi": (data_dict[case]["VM_S12"][self.num][0] + data_dict[case]["VM_S12"][self.num][1]) / 2,
                            "S12": (data_dict[case]["VM_S12"][self.num][2] + data_dict[case]["VM_S12"][self.num][3]) / 2
                            }
                     for case in data_dict
                     }

        self.data_norm = {case: {"VMi": (self.data[case]["VMi"] - VMi_mm[case][0]) / (VMi_mm[case][1] - VMi_mm[case][0]),
                                 "S12": (self.data[case]["S12"] - S12_mm[case][0]) / (S12_mm[case][1] - S12_mm[case][0])}
                          for case in self.data
                          }

    def __repr__(self):
        return f'Element {self.num}: nodes: \n{self.corners}\n'

    def __lt__(self, other):
        if isinstance(other, Element):
            return min(c.x for c in self.corners) < min(c.x for c in other.corners)
        else:
            raise TypeError(f"Cannot compare Element to {type(other)}")

    def plot_self(self, ax, case, tpe):
        self.xs, self.ys, self.zs = np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))
        for i, corner in enumerate(self.corners):
            corner.displace(case)
            self.xs[i % 2, i // 2], self.ys[i % 2, i // 2], self.zs[i % 2, i // 2] = corner.xd, corner.yd, corner.zd

        cmp = plt.get_cmap('viridis')
        col = [[cmp(self.data_norm[case][tpe])]]

        ax.plot_surface(self.xs, self.zs, self.ys, facecolors=col)

    def find_2d(self, case):
        z, y, x = [], [], []
        for i, corner in enumerate(self.corners):
            corner.displace(case)
            if corner.zd not in z and corner.yd not in y:
                z.append(corner.zd)
                y.append(corner.yd)
            if corner.x not in x:
                x.append(corner.x)

        self.z, self.y, self.x = z, y, x

    def plot_2d(self, ax, case, tpe):
        self.find_2d(case)

        cmp = plt.get_cmap('viridis')
        col = cmp(self.data_norm[case][tpe])

        ax.plot(self.z, self.y, color=col)


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
        return f'Node {self.num}: ({self.x}, {self.y}, {self.z})\n'

    def __lt__(self, other):
        if isinstance(other, Node):
            return self.x < other.x
        else:
            raise TypeError(f"Cannot compare Node to {type(other)}")

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


def read_inp(data_dict, VMi_mm, S12_mm):
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
        node_lst = []

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

                elem_lst.append(Element(data_line[0], data_line[1:], node_dct, data_dict, VMi_mm, S12_mm))

            if nodes and s:
                data_line = line.strip("\n").replace(' ', '')
                data_line = data_line.split(",")
                for i, data_pt in enumerate(data_line):
                    data_line[i] = float(data_pt) if i > 0 else int(data_pt)

                node_inst = Node(data_line[0], tuple(data_line[1:]), data_dict, False)
                node_dct[data_line[0]] = node_inst
                node_lst.append(node_inst)

            s = True

    return node_dct, node_lst, elem_lst, asem_lst


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
            VMi_lst[case] = [(data_dict[case]["VM_S12"][elm][0] + data_dict[case]["VM_S12"][elm][1]) / 2 for elm in data_dict[case]["VM_S12"]]
            S12_lst[case] = [(data_dict[case]["VM_S12"][elm][2] + data_dict[case]["VM_S12"][elm][3]) / 2 for elm in data_dict[case]["VM_S12"]]

        VMi = {case: (min(VMi_lst[case]), max(VMi_lst[case])) for case in VMi_lst}
        S12 = {case: (min(S12_lst[case]), max(S12_lst[case])) for case in S12_lst}

    return data_dict, VMi, S12


if __name__ == '__main__':
    ipt = input("Load case and stress to display: ")
    case, tpe = ipt.split(", ")
    data, VMi, S12 = read_rpt()
    node_dict, node, elem, asem = read_inp(data, VMi, S12)

    mm = VMi[case] if tpe == "VMi" else S12[case]

    fig = plt.figure()
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)
    ax1.set_ylim(0 - 600, 100)
    ax1.set_zlim(0 - 350, 350)
    ax2.set_xlim(0 - 550, 150)
    ax2.set_ylim(0 - 350, 350)

    node.sort()
    elem.sort()

    [el_min] = [(i, el) for i, el in enumerate(elem) if el.data[case][tpe] == mm[0]]
    [el_max] = [(i, el) for i, el in enumerate(elem) if el.data[case][tpe] == mm[1]]
    i = (el_min[0] // 62) * 62
    j = (el_max[0] // 62) * 62

    for el in elem[j:j+62]:
        el.plot_2d(ax2, case, tpe)

    x_mdl = round(sum(el.x) / len(el.x), 2)


    for el in elem:
        el.plot_self(ax1, case, tpe)

    for no in asem:
        scale = 100 if no.num in range(5, 16) else 1
        no.plot_quiver(ax1, case, u=True, sc=scale)

    normal = color.Normalize(vmin=mm[0], vmax=mm[1])
    cmap = cm.ScalarMappable(norm=normal, cmap=plt.viridis())
    cbar = plt.colorbar(cmap)

    cbl = "Von Mises stress" if tpe == "VMi" else "Shear stress"
    cbar.set_label(cbl)

    if case == "Bending":
        ttl = "Bending case"
    elif case == "Jam_Bent":
        ttl = "Bending and jammed actuator case"
    elif case == "Jam_Straight":
        ttl = "Only jammed actuator case"
    else:
        ttl = ""
    fig.suptitle(f"{cbl} distribution in the deformed aileron\n{ttl}")
    ax1.set_title("Entire Aileron")
    ax2.set_title(f"Cross section with maximum stress at x={x_mdl}mm")

    ax1.set_xlabel("X")
    ax1.set_ylabel("Z")
    ax1.set_zlabel("Y")
    ax2.set_xlabel("Z")
    ax2.set_ylabel("Y")
    plt.show()
