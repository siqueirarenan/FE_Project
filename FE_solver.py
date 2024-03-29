import numpy as np
import FE_classes
from scipy import linalg as spln
import pickle


# Stiffness matrix construction

def solver(job_name, mdb):
    mdl = mdb.models[mdb.jobs[job_name].model]
    for p in mdl.parts.values():  # Construction of u general vector
        # Stiffness Matrix construction
        elementsBmatrix = {}
        elementsEmatrix = {}
        elementsstiffness = {}
        for e in p.elements:  # Begining of construction of K_e
            elementsEmatrix[e.label] = ElementMaterialMatrix(mdl, e)
            elementsBmatrix[e.label], detJ = ElementBmatrix(p, e)
            elementsstiffness[e.label] = ElementStiffness(elementsBmatrix[e.label], elementsEmatrix[e.label], detJ)
        # Assembly
        K = np.zeros((3 * len(p.nodes), 3 * len(p.nodes)))
        for k_e in elementsstiffness.items():
            for m in range(p.elements.nodePerElement):
                for n in range(p.elements.nodePerElement):  # Only up half of the simmetrtic matrix is calculated
                    K[3 * p.elements[k_e[0] - 1].connectivity[m]:3 * p.elements[k_e[0] - 1].connectivity[m] + 3,
                    3 * p.elements[k_e[0] - 1].connectivity[n]:3 * p.elements[k_e[0] - 1].connectivity[n] + 3
                    ] += k_e[1][3 * m:3 * m + 3, 3 * n:3 * n + 3]
        # Step
        for s in mdl.steps.values():  # per Step
            # Load vector
            F = np.zeros((3 * len(p.nodes), 1))
            for f in s.loadStates.items():  # Checking each load in the step
                if isinstance(f[1], FE_classes.ConcentratedForceState):  # Dealing according to case
                    gen_force = np.array([f[1].cf1, f[1].cf2, f[1].cf3]) / len(
                        mdl.loads[f[0]].region.nodes)
                for n in mdl.loads[f[0]].region.nodes:  # checking each node of the load
                    F[3 * (n.label - 1):3 * (n.label - 1) + 3, 0] += gen_force

            # Boundary conditions
            # u = np.empty((3 * len(p.nodes), 1))
            rows_columns_toNull = []
            for bc in s.boundaryConditionStates.items():  # Checking each BC in the step
                if isinstance(bc[1], FE_classes.DisplacementBCState):  # Dealing according to case
                    gen_bc = np.array([bc[1].u1, bc[1].u2, bc[1].u3])
                for n in mdl.boundaryConditions[bc[0]].region.nodes:  # checking each node of the load
                    # u[3 * (n.label-1):3 * (n.label-1) + 3 ,0] +=  gen_bc #still not implemented for values != 0
                    # getting the doF with null displacement
                    for ii in range(3):
                        rows_columns_toNull += [3 * (n.label - 1) + ii] if gen_bc[ii] == 0 else []
            # Eliminating rows and columns
            K = np.delete(np.delete(K, rows_columns_toNull, 0), rows_columns_toNull, 1)
            F = np.delete(F, rows_columns_toNull, 0)
            # SOLVER Ku=F
            u = spln.solve(K, F, assume_a='sym')
            # Assembling u again
            uu = np.zeros((3 * len(p.nodes)))
            count = 0
            for dof in range(3 * len(p.nodes)):
                if any([dof == rct for rct in rows_columns_toNull]):
                    uu[dof] = 0
                else:
                    uu[dof] = u[count]
                    count += 1
            uu = uu.reshape((len(p.nodes), 3))

            ## ODB creation
            odb = FE_classes.Odb('Output_files\\' + job_name)
            initial_frame = odb.Step(s.name.upper()).Frame()
            frame = odb.Step(s.name.upper()).Frame()
            _ = frame.FieldOutput('U').values.NodeArray(uu)
            instance = odb.rootAssembly.Instance(p.name.upper() + '-1', p)
            for sets in instance.obj.sets.values():  # Distributing Mdb sets into Odb sets
                if len(sets.elements) > 0:
                    instance.ElementSet(sets.name.upper(), sets.elements)
            historyregion = odb.steps[s.name.upper()].HistoryRegion('Assembly ASSEMBLY')

            # Field outputs
            for fo in s.fieldOutputRequestStates.values():
                for v in fo.variables:
                    if v == 'MISESMAX':
                        vonmises = []
                        for e in p.elements:
                            u_e = []
                            for c in e.connectivity:
                                u_e = np.hstack([u_e, uu[c]])
                            strain_e = elementsBmatrix[e.label] @ u_e
                            s_e = elementsEmatrix[e.label] @ strain_e
                            vonmises += [(((s_e[0] - s_e[1]) ** 2 + (s_e[1] - s_e[2]) ** 2 + (s_e[2] -
                                                                                              s_e[0]) ** 2 + 6 * (
                                                       s_e[3] ** 2 + s_e[4] ** 2 + s_e[5] ** 2)) / 2) ** 0.5]
                            frame.FieldOutput('MISESMAX').values.ElementArray(vonmises)
                    elif v == 'ESEDEN':
                        eseden = []
                        for e in p.elements:
                            u_e = []
                            for c in e.connectivity:
                                u_e = np.hstack([u_e, uu[c]])
                            eseden += [(u_e.T @ elementsstiffness[e.label] @ u_e) / 2]
                            frame.FieldOutput('ESEDEN').values.ElementArray(eseden)
                    elif v == 'EVOL':
                        if p.elementType == 'uniformHex':
                            evol = [p.HexElementVolume for e in p.elements]
                        frame.FieldOutput('EVOL').values.ElementArray(evol)
                    else:
                        print('Field output {} not found'.format(v))
            # History outputs
            for ho in s.historyOutputRequestStates.values():
                for v in ho.variables:
                    if v == 'ALLWK':  # Total energy output
                        histout = historyregion.HistoryOutput(v, initial_frame.number, 0)
                        compliance_energy = F.T @ u
                        histout.addData(frame.number, compliance_energy[0, 0])
                    else:
                        print('History output {} not found'.format(v))

    with open('Output_files\\' + job_name + '.odb', 'wb') as output:
        pickle.dump(odb, output, pickle.HIGHEST_PROTOCOL)


def myInverse(A):
    detA = A[0, 0] * A[1, 1] * A[2, 2] + A[1, 0] * A[2, 1] * A[0, 2] + A[2, 0] * A[0, 1] * A[1, 2] - A[0,
                                                                                                       0] * A[2, 1] * A[
               1, 2] - A[2, 0] * A[1, 1] * A[0, 2] - A[1, 0] * A[0, 1] * A[2, 2]
    b00 = A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]
    b01 = A[0, 2] * A[2, 1] - A[0, 1] * A[2, 2]
    b02 = A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]
    b10 = A[1, 2] * A[2, 0] - A[1, 0] * A[2, 2]
    b11 = A[0, 0] * A[2, 2] - A[0, 2] * A[2, 0]
    b12 = A[0, 2] * A[1, 0] - A[0, 0] * A[1, 2]
    b20 = A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]
    b21 = A[0, 1] * A[2, 0] - A[0, 0] * A[2, 1]
    b22 = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    Ainv = np.array([[b00, b01, b02], [b10, b11, b12], [b20, b21, b22]]) / detA

    return Ainv, detA


def ElementBmatrix(partObj, elementObj):
    u_e = np.array(
        partObj.nodes[elementObj.connectivity[0]].coordinates)  # Construction of u_e = [x1 y1 z1; x2 y2 z2; x3... ]
    for n in elementObj.connectivity[1:]:
        u_e = np.vstack([u_e, partObj.nodes[n].coordinates])
    # Constant DQ - derivative of the 8 N functions from i, j and k (rows); after replacing i,j,k = 0
    DQ = 0.125 * np.array([[-1., 1., -1., 1., -1., 1., -1., 1.],
                           [-1., -1., 1., 1., -1., -1., 1., 1.],
                           [-1., -1., -1., -1., 1., 1., 1., 1.]])
    # Jacobian inverse already subs zero
    J_inv, detJ = myInverse(DQ.dot(u_e))
    # Matrics constining the derivatives of N through x,y,z, already with i,j,s = 0
    dQ = np.zeros((partObj.elements.nodePerElement, 3))
    for n in range(partObj.elements.nodePerElement):
        dQ[n, :] = np.array(np.matrix.transpose(J_inv.dot(DQ[:, n])))

    # Strin-displacemente matrix
    B = np.array([[dQ[0, 0], 0, 0, dQ[1, 0], 0, 0, dQ[2, 0], 0, 0, dQ[3, 0], 0, 0, dQ[4, 0], 0, 0, dQ[5, 0], 0, 0,
                   dQ[6, 0], 0, 0, dQ[7, 0], 0, 0],
                  [0, dQ[0, 1], 0, 0, dQ[1, 1], 0, 0, dQ[2, 1], 0, 0, dQ[3, 1], 0, 0, dQ[4, 1], 0, 0, dQ[5, 1], 0, 0,
                   dQ[6, 1], 0, 0, dQ[7, 1], 0],
                  [0, 0, dQ[0, 2], 0, 0, dQ[1, 2], 0, 0, dQ[2, 2], 0, 0, dQ[3, 2], 0, 0, dQ[4, 2], 0, 0, dQ[5, 2], 0, 0,
                   dQ[6, 2], 0, 0, dQ[7, 2]],
                  [0, dQ[0, 2], dQ[0, 1], 0, dQ[1, 2], dQ[1, 1], 0, dQ[2, 2], dQ[2, 1], 0, dQ[3, 2], dQ[3, 1], 0,
                   dQ[4, 2], dQ[4, 1], 0, dQ[5, 2], dQ[5, 1], 0, dQ[6, 2], dQ[6, 1], 0, dQ[7, 2], dQ[7, 1]],
                  [dQ[0, 2], 0, dQ[0, 0], dQ[1, 2], 0, dQ[1, 0], dQ[2, 2], 0, dQ[2, 0], dQ[3, 2], 0, dQ[3, 0], dQ[4, 2],
                   0, dQ[4, 0], dQ[5, 2], 0, dQ[5, 0], dQ[6, 2], 0, dQ[6, 0], dQ[7, 2], 0, dQ[7, 0]],
                  [dQ[0, 1], dQ[0, 0], 0, dQ[1, 1], dQ[1, 0], 0, dQ[2, 1], dQ[2, 0], 0, dQ[3, 1], dQ[3, 0], 0, dQ[4, 1],
                   dQ[4, 0], 0, dQ[5, 1], dQ[5, 0], 0, dQ[6, 1], dQ[6, 0], 0, dQ[7, 1], dQ[7, 0], 0]])
    return B, detJ


def ElementMaterialMatrix(mdlObj, elementObj):
    # Element material
    E = mdlObj.materials[elementObj.material].elastic.table[0]
    v = mdlObj.materials[elementObj.material].elastic.table[1]
    p1 = v * E / ((1 + v) * (1 - 2 * v))
    p2 = E / (2 * (1 + v))
    E_e = np.array([[p1 + 2 * p2, p1, p1, 0, 0, 0],
                    [p1, p1 + 2 * p2, p1, 0, 0, 0],
                    [p1, p1, p1 + 2 * p2, 0, 0, 0],
                    [0, 0, 0, p2, 0, 0],
                    [0, 0, 0, 0, p2, 0],
                    [0, 0, 0, 0, 0, p2]])
    return E_e


def ElementStiffness(B, E_e, detJ):
    K_e = 2 * 2 * 2 * detJ * B.T @ E_e @ B
    correc = np.array([[9.93189967e-01, 1.00000000e+00, 1.00000000e+00, 9.87585393e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00344016e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 9.93189967e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 9.87585393e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00344016e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 9.93189967e-01, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 9.87585393e-01, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00344016e+00],
                       [9.87585393e-01, 1.00000000e+00, 1.00000000e+00, 9.93189967e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00344016e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        9.93189967e-01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 9.87585393e-01, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00344016e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 9.93189967e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 9.87585393e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00344016e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 9.93189967e-01, 1.00000000e+00,
                        1.00000000e+00, 9.87585393e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00344016e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 9.87585393e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 9.93189967e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00344016e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        9.93189967e-01, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00344016e+00, 1.00000000e+00, 1.00000000e+00,
                        9.87585393e-01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 9.87585393e-01, 1.00000000e+00,
                        1.00000000e+00, 9.93189967e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00344016e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        9.87585393e-01, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 9.93189967e-01, 1.00000000e+00,
                        1.00000000e+00, 1.00344016e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 9.93189967e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00344016e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 9.87585393e-01],
                       [1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00344016e+00, 1.00000000e+00, 1.00000000e+00,
                        9.93189967e-01, 1.00000000e+00, 1.00000000e+00, 9.87585393e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00344016e+00, 1.00000000e+00,
                        1.00000000e+00, 9.93189967e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 9.87585393e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 9.87585393e-01, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00344016e+00,
                        1.00000000e+00, 1.00000000e+00, 9.93189967e-01, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00344016e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        9.87585393e-01, 1.00000000e+00, 1.00000000e+00, 9.93189967e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00344016e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        9.93189967e-01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 9.87585393e-01, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 9.87585393e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00344016e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 9.93189967e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00344016e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 9.93189967e-01, 1.00000000e+00,
                        1.00000000e+00, 9.87585393e-01, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00344016e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 9.87585393e-01, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 9.93189967e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00344016e+00, 1.00000000e+00, 1.00000000e+00,
                        9.87585393e-01, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        9.93189967e-01, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00],
                       [1.00344016e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 9.87585393e-01, 1.00000000e+00,
                        1.00000000e+00, 9.93189967e-01, 1.00000000e+00, 1.00000000e+00],
                       [1.00000000e+00, 1.00344016e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00541661e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        9.87585393e-01, 1.00000000e+00, 1.00000000e+00, 1.00541661e+00,
                        1.00000000e+00, 1.00000000e+00, 9.93189967e-01, 1.00000000e+00],
                       [1.00000000e+00, 1.00000000e+00, 1.00344016e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 9.87585393e-01,
                        1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00000000e+00, 1.00541661e+00, 1.00000000e+00, 1.00000000e+00,
                        1.00541661e+00, 1.00000000e+00, 1.00000000e+00, 9.93189967e-01]])
    K_e = np.divide(K_e, correc)

    return K_e
