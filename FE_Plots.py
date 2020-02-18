from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
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
def deformedNodePlot(mdl, p, s, uu, scale_factor=1):
    fig = plt.figure()
    uu = uu.reshape((len(p.nodes), 3))
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

def vonMisesHexMeshPlot(w,h,d,vonmises):
    from matplotlib.colors import to_hex
    colornew = []
    for i in range(256):
        c = 2 * i / 256 - 1
        if c < -0.8:
            b = (((c + 1) / (-0.8 + 1)) * (1 - 0.6)) + 0.6
            r, g = 0, 0
        elif c < -0.25:
            g = (((c + 0.8) / (-0.25 + 0.8)) * (1))
            b, r = 1, 0
        elif c < 0.25:
            g = 1
            r = (((c + 0.25) / (0.25 + 0.25)) * (1))
            b = (((c + 0.25) / (0.25 + 0.25)) * (-1)) + 1
        elif c < 0.8:
            g = (((c - 0.25) / (0.8 - 0.25)) * (-1)) + 1
            b, r = 0, 1
        else:
            r = (((c - 0.8) / (1 - 0.8)) * (0.8 - 1)) + 1
            b, g = 0, 0
        colornew += [[r, g, b]]
    vm = np.array([list(vonmises.values())])
    vm = np.divide(vm,max(max(vm)))
    vm = vm.reshape((d,h,w))
    color = np.empty((d,h,w), dtype = 'object')
    for x in range(vm.shape[0]):
        for y in range(vm.shape[1]):
            for z in range(vm.shape[2]):
                color[x,y,z] = to_hex(colornew[int(round(vm[x,y,z]*len(colornew))-1)])
    fig = plt.figure()
    ax = fig.add_subplot('111', projection='3d')
    ax.voxels(np.ones([d, h, w], dtype=np.bool), facecolors=color, edgecolors='k')