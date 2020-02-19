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
    def get_spar(x, y, z, n, sort=False):
        c = np.where(z == 0)
        xsp, ysp, zsp, nsp = x[c], y[c], z[c], n[c]

        spar = {}
        for i, nd in enumerate(nsp):
            spar[nd] = (xsp[i], ysp[i], zsp[i])

        if sort:
            zsort = np.argsort(zsp)
            xsp, ysp, zsp, nsp = xsp[zsort], ysp[zsort], zsp[zsort], nsp[zsort]
            xsort = np.argsort(xsp)
            xsp, ysp, zsp, nsp = xsp[xsort], ysp[xsort], zsp[xsort], nsp[xsort]

        return xsp, ysp, zsp, nsp, spar

    def get_top(x, y, z, n, sort=False):
        c = np.where(y >= -1.05331139e-08)
        xtp, ytp, ztp, ntp = x[c], y[c], z[c], n[c]

        _, ysp, _, nsp, _ = get_spar(x, y, z, n)
        nsp = nsp[np.where(ysp != max(ysp))]
        c = np.where(np.isin(ntp, nsp))
        xtp, ytp, ztp, ntp = np.delete(xtp, c), np.delete(ytp, c), np.delete(ztp, c), np.delete(ntp, c)

        top = {}
        for i, nd in enumerate(ntp):
            top[nd] = (xtp[i], ytp[i], ztp[i])

        if sort:
            xsort = np.argsort(xtp)
            xtp, ytp = xtp[xsort], ytp[xsort]
            zsort = np.argsort(ztp)
            ztp, ytp = ztp[zsort], ytp[zsort]

        return xtp, ytp, ztp, ntp, top

    def get_bottom(x, y, z, n):
        c = np.where(y <= 0)
        xbp, ybp, zbp, nbp = x[c], y[c], z[c], n[c]

        _, ysp, _, nsp, _ = get_spar(x, y, z, n)
        nsp = nsp[np.where(ysp != min(ysp))]
        c = np.where(np.isin(nbp, nsp))
        xbp, ybp, zbp, nbp = np.delete(xbp, c), np.delete(ybp, c), np.delete(zbp, c), np.delete(nbp, c)

        bottom = {}
        for i, nd in enumerate(nbp):
            bottom[nd] = (xbp[i], ybp[i], zbp[i])

        return xbp, ybp, zbp, nbp, bottom

    def get_nodes(nodes):
        x = np.zeros((len(nodes),))
        y = np.zeros((len(nodes),))
        z = np.zeros((len(nodes),))
        n = np.zeros((len(nodes),))

        for node_num, node_coord in node_dict.items():
            x[node_num - 1] = node_coord[0]
            y[node_num - 1] = node_coord[1]
            z[node_num - 1] = node_coord[2]
            n[node_num - 1] = node_num

        return x, y, z, n, nodes

    def find_shape_xz(x, y, z, n):
        xco = []
        zco = []

        for i in range(len(n)):
            if x[i] not in xco:
                xco.append(x[i])
            if z[i] not in zco:
                zco.append(z[i])

        return len(zco), len(xco)

    fig = plt.figure()
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    xn, yn, zn, nn, _ = get_nodes(node_dict)
    xs, ys, zs, ns, _ = get_spar(xn, yn, zn, nn, sort=True)
    xt, yt, zt, nt, _ = get_top(xn, yn, zn, nn, sort=True)
    xb, yb, zb, nb, _ = get_bottom(xn, yn, zn, nn)

    #print(len(yt[np.where(xt == 172)]))
    #print(yt[np.where(xt == 172)])
    #print(sorted(zt[np.where(xt == 172)]))
    shapet = find_shape_xz(xt, yt, zt, nt)
    print(shapet)
    np.reshape(yt, shapet)

    img = ax1.scatter(xt, zt, yt, c=nt, cmap=plt.hot())
    ax2.plot_wireframe(np.reshape(xt, shapet), np.reshape(zt, shapet), np.reshape(yt, shapet))
    ax1.set_xlabel("X")
    ax1.set_ylabel("Z")
    ax1.set_zlabel("Y")
    fig.colorbar(img)

    plt.show()


if __name__ == '__main__':
    node = read_inp()
    plot_aileron(node)
