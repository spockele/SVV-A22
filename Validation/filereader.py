"""
Reader for the .rpt and .inp files
"""

import matplotlib.pyplot as plt
import matplotlib.colors as color
import matplotlib.cm as cm
import mpl_toolkits.mplot3d as axes3d
import numpy as np
import tabulate as tab


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

    def plot_U(self, ax, case):
        if case == "Bending":
            _, _, dy, dz = self.U_Bend
        elif case == "Jam_Bent":
            _, _, dy, dz = self.U_JBen
        elif case == "Jam_Straight":
            _, _, dy, dz = self.U_JStr

        ax.plot([self.x/1000], [dy/1000], "ro")
        ax.plot([self.x/1000], [dz/1000], "bo")


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


def table_reaction(data_dict, case):
    def extract_relevant():
        data_rel = []
        for i in range(3):
            h = (1, 3, 16)[i]
            data_rel.append([f"Hinge {i+1}"] + [round(n, 3) for n in data_dict[case]["reaction"][h]])

        data_rel.append([f"Actuator 1"] + [round(n, 3) for n in data_dict[case]["reaction"][2]])

        return data_rel

    def save_table(tab):
        with open(f"Table_{case}.txt", "w") as f:
            f.writelines(tab)

    data_tab = extract_relevant()
    table = tab.tabulate(data_tab, ("Point", "|F|", "Fx", "Fy", "Fz"), numalign="right", tablefmt="presto")
    latex = tab.tabulate(data_tab, ("Point", "|F|", "Fx", "Fy", "Fz"), numalign="right", tablefmt="latex_raw")
    save_table(latex)
    print(table)


def main(case, tpe, view):
    # Pre-Processing of the data files
    print(f"\nPlot_{case}_{tpe}_{view}\n-----------------------------------------------------")
    data, VMi, S12 = read_rpt()
    node_dict, node, elem, asem = read_inp(data, VMi, S12)

    table_reaction(data, case)

    mm = VMi[case] if tpe == "VMi" else S12[case]

    node.sort()
    elem.sort()

    [el_min, *_] = [(i, el) for i, el in enumerate(elem) if el.data[case][tpe] == mm[0]]
    [el_max, *_] = [(i, el) for i, el in enumerate(elem) if el.data[case][tpe] == mm[1]]
    idx = (el_min[0] // 62) * 62 if abs(el_min[1].data[case][tpe]) > abs(el_max[1].data[case][tpe])\
        else (el_max[0] // 62) * 62

    # Plotting of the data
    fig = plt.figure()
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122)
    ax1.set_ylim(0 - 600, 100)
    ax1.set_zlim(0 - 350, 350)
    ax2.set_xlim(0 - 550, 150)
    ax2.set_ylim(0 - 350, 350)

    for el in elem[idx:idx+62]:
        el.plot_2d(ax2, case, tpe)

    x_mdl = round(sum(elem[idx].x) / len(elem[idx].x), 2)

    for el in elem:
        el.plot_self(ax1, case, tpe)

    # for no in asem:
    #     scale = 100 if no.num in range(5, 16) else 1
    #     no.plot_quiver(ax1, case, u=True, sc=scale)

    normal = color.Normalize(vmin=mm[0], vmax=mm[1])
    cmap = cm.ScalarMappable(norm=normal, cmap=plt.viridis())
    cbar = plt.colorbar(cmap)

    cbt = "Von Mises stress [GPa]" if tpe == "VMi" else "Shear stress [GPa]"
    cbar.set_label(cbt)

    if case == "Bending":
        ttl = "Bending case"
    elif case == "Jam_Bent":
        ttl = "Bending and jammed actuator case"
    elif case == "Jam_Straight":
        ttl = "Only jammed actuator case"
    else:
        ttl = ""

    ax1.set_xlabel("X [mm]")
    ax1.set_ylabel("Z [mm]")
    ax1.set_zlabel("Y [mm]")
    ax2.set_xlabel("Z [mm]")
    ax2.set_ylabel("Y [mm]")

    plt.subplots_adjust(0.05, 0.05, 0.95, 0.95, 0.15, 0.15)
    if view == "top":
        ax1.view_init(azim=40, elev=30)
    elif view == "bottom":
        ax1.view_init(azim=40, elev=-30)
    else:
        ax1.view_init(azim=40, elev=0)

    plt.show()

    le = (1,   2,   5,   7,   9,  11,  12,  43,  44,  57,  58,  59,  60,  61,  62,  63,
          64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
          80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
          96,  97,  98,  99, 100, 101, 102, 103, 156, 157, 158, 159, 160, 161, 173, 174,
          175, 176, 177, 178, 190, 191, 192, 193, 194, 195, 212, 213, 214, 215, 216, 217,
          218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
          234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245
          )

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    [no.plot_U(ax1, case) for n, no in node_dict.items() if n in le]
    ax1.set_xlabel("x [m]")
    ax1.set_ylabel("Deflection [m]")
    ax1.legend(["Deflection in y", "Deflection in x"])
    plt.show()


if __name__ == '__main__':
    main("Jam_Bent", "VMi", "top")
    main("Jam_Bent", "S12", "top")
    #for case in ("Bending", "Jam_Bent", "Jam_Straight"):
    #    for tpe in ("VMi", "S12"):
    #        for view in ("top", "bottom"):
    #            main(case, tpe, view)
