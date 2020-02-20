from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
from matplotlib.collections import EventCollection
import numpy as np


# Undeformed nodes showing loads and BC
def undeformedNodePlot(mdl, p, s):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for n in p.nodes:  # Ploting nodes
        ax.scatter(n.coordinates[0], n.coordinates[1],
                   n.coordinates[2], marker='.', color=[0.3, 0.3, 0.3])
    for l in s.loadStates.items():  # Marking load nodes
        for n in mdl.loads[l[0]].region.nodes:
            ax.scatter(n.coordinates[0], n.coordinates[1],
                       n.coordinates[2], marker='o', color='blue')
    for bc in s.boundaryConditionStates.items():  # Marking BC nodes
        for n in mdl.boundaryConditions[bc[0]].region.nodes:
            ax.scatter(n.coordinates[0], n.coordinates[1],
                       n.coordinates[2], marker='o', color='red')


# Deformed nodes
def deformedNodePlot(mdl, p, s, odb, scale_factor=1):
    fig = plt.figure()
    uu = odb.steps[s.name].frames[-1].fieldOutputs['U'].values.original_data
    ug = scale_factor * uu
    ax = fig.add_subplot(111, projection='3d')
    n_size_o = []
    for n in p.nodes:  # Calculate the norm of the displacements
        n_size_o += [np.linalg.norm(uu[n.label - 1, :])]
    n_size = n_size_o / max(n_size_o)
    for n in p.nodes:
        ax.scatter(n.coordinates[0] + ug[n.label - 1, 0], n.coordinates[1] + ug[n.label - 1, 1],
                   n.coordinates[2] + ug[n.label - 1, 2], marker='.',
                   color=np.array([n_size[n.label - 1], 0, 1 - n_size[n.label - 1]]))
    for bc in s.boundaryConditionStates.items():  # Marking BC nodes
        for n in mdl.boundaryConditions[bc[0]].region.nodes:
            ax.scatter(n.coordinates[0] + ug[n.label - 1, 0], n.coordinates[1] + ug[n.label - 1, 1],
                       n.coordinates[2] + ug[n.label - 1, 2], marker='o', color='red')


def FieldOutputHexMeshPlot(p, s, odb, fieldoutputname, sets=None):
    from matplotlib.colors import to_hex
    w, h, d = p.HexElementSizes
    x, y, z = w, h, d
    fieldoutput = odb.steps[s.name.upper()].frames[-1].fieldOutputs[fieldoutputname].values.original_data
    colornew = []
    for i in range(256):
        c = 2 * i / 256 - 1
        if c < -0.8:
            b = (((c + 1) / (-0.8 + 1)) * (1 - 0.6)) + 0.6
            r, g = 0, 0
        elif c < -0.25:
            g = (((c + 0.8) / (-0.25 + 0.8)) * 1)
            b, r = 1, 0
        elif c < 0.25:
            g = 1
            r = (((c + 0.25) / (0.25 + 0.25)) * 1)
            b = (((c + 0.25) / (0.25 + 0.25)) * (-1)) + 1
        elif c < 0.8:
            g = (((c - 0.25) / (0.8 - 0.25)) * (-1)) + 1
            b, r = 0, 1
        else:
            r = (((c - 0.8) / (1 - 0.8)) * (0.8 - 1)) + 1
            b, g = 0, 0
        colornew += [[r, g, b]]
    vm = np.array(fieldoutput)
    vm = np.divide((vm - min(vm)), max(vm - min(vm)))
    vm = vm.reshape((d, h, w))
    color = np.empty(p.HexElementSizes, dtype='object')
    for xx in range(vm.shape[0]):
        for yy in range(vm.shape[1]):
            for zz in range(vm.shape[2]):
                color[xx, yy, zz] = to_hex(colornew[int(round(vm[zz, yy, xx] * (len(colornew) - 1)))])
    if sets is None:
        set_regions = np.ones(p.HexElementSizes, dtype=np.bool)
    else:
        set_regions = np.zeros(p.HexElementSizes, dtype=bool)
        for s_set in sets:
            set_region = np.zeros(p.HexElementSizes, dtype=bool)
            for e in p.sets[s_set].elements:
                xe, ye, ze = e.label % x, (e.label // x + 1) % y, ((e.label // x + 1) // y + 1) % z
                set_region[xe - 1, ye - 1, ze - 1] = True
            set_regions |= set_region
    fig = plt.figure()
    ax = fig.add_subplot('111', projection='3d')
    ax.voxels(set_regions, facecolors=color, edgecolors='k')
    ax.set_xlim3d(int((x - max(x, y, z)) / 2) - 1, int((x + max(x, y, z)) / 2) + 1)
    ax.set_ylim3d(int((y - max(x, y, z)) / 2) - 1, int((y + max(x, y, z)) / 2) + 1)
    ax.set_zlim3d(int((z - max(x, y, z)) / 2) - 1, int((z + max(x, y, z)) / 2) + 1)
    plt.axis('off')
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)


def HistoryOutputPlot(odb, s, ho_names):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for ho in ho_names:
        data = odb.steps[s.name.upper()].historyRegions['Assembly ASSEMBLY'].historyOutputs[ho].data
        xdata, ydata = [], []
        for d in data:
            xdata += [d[0]]
            ydata += [d[1]]
        ax.plot(xdata, ydata, color='tab:blue')


def HexGeometryPlot(fig, p, sets):
    x, y, z = p.HexElementSizes
    set_regions = np.zeros(p.HexElementSizes, dtype=bool)
    colors = np.empty(p.HexElementSizes, dtype=object)
    color_sequence = ('grey', 'blue', 'green', 'red', 'yellow',)
    count = 0
    for s_set in sets:
        set_region = np.zeros(p.HexElementSizes, dtype=bool)
        for e in p.sets[s_set].elements:
            xe, ye, ze = e.label % x, (e.label // x + 1) % y, ((e.label // x + 1) // y + 1) % z
            set_region[xe - 1, ye - 1, ze - 1] = True
        set_regions |= set_region
        colors[set_region] = color_sequence[count]
        count += 1
    ax = fig.add_subplot('111', projection='3d')
    ax.voxels(set_regions, facecolors=colors, edgecolor='k')
    ax.set_xlim3d(int((x - max(x, y, z)) / 2) - 1, int((x + max(x, y, z)) / 2) + 1)
    ax.set_ylim3d(int((y - max(x, y, z)) / 2) - 1, int((y + max(x, y, z)) / 2) + 1)
    ax.set_zlim3d(int((z - max(x, y, z)) / 2) - 1, int((z + max(x, y, z)) / 2) + 1)
    plt.axis('off')
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
