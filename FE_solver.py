import numpy as np
import FE_classes
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import sympy as sp

#Stiffness matrix construction

def solver(job_name,mdb):

    for p in mdb.models[mdb.jobs[job_name].model].parts.values(): #Construction of u general vector
        #Stiffness Matrix construction
        elementsstiffness = {}
        for e in p.elements:                                 #Begining of construction of K_e
            elementsstiffness[e.label] = ElementStiffness(mdb.models[mdb.jobs[job_name].model],p,e)
        #Assembly
        K = np.empty((3*len(p.nodes),3*len(p.nodes)))
        for k_e in elementsstiffness.items():
            for m in range(p.elements.nodePerElement):
                for n in range(p.elements.nodePerElement): #Only up half of the simmetrtic matrix is calculated
                    K[3 * p.elements[k_e[0]-1].connectivity[m]:3 * p.elements[k_e[0]-1].connectivity[m] + 3,
                      3 * p.elements[k_e[0]-1].connectivity[n]:3 * p.elements[k_e[0]-1].connectivity[n] + 3
                    ] += k_e[1][3 * m:3*m + 3, 3 * n:3*n + 3]
        #Step
        for s in mdb.models[mdb.jobs[job_name].model].steps.values(): #per Step
            #Load vector
            F = np.empty((3 * len(p.nodes),1))
            for f in s.loadStates.items(): #Checking each load in the step
                if isinstance(f[1],FE_classes.ConcentratedForceState): #Dealing according to case
                    gen_force = np.array([f[1].cf1,f[1].cf2 ,f[1].cf3]) / len(
                        mdb.models[mdb.jobs[job_name].model].loads[f[0]].region.nodes)
                for n in mdb.models[mdb.jobs[job_name].model].loads[f[0]].region.nodes:  #checking each node of the load
                    F[3 * (n.label-1):3 * (n.label-1) + 3 ,0] +=  gen_force
            #Boundary conditions
            #u = np.empty((3 * len(p.nodes), 1))
            rows_columns_toNull = []
            for bc in s.boundaryConditionStates.items():  # Checking each BC in the step
                if isinstance(bc[1], FE_classes.DisplacementBCState):  # Dealing according to case
                    gen_bc = np.array([bc[1].u1,bc[1].u2 ,bc[1].u3])
                for n in mdb.models[mdb.jobs[job_name].model].boundaryConditions[bc[0]].region.nodes:  #checking each node of the load
                    #u[3 * (n.label-1):3 * (n.label-1) + 3 ,0] +=  gen_bc #still not implemented for values != 0
                    #getting the doF with null displacement
                    for ii in range(3):
                        rows_columns_toNull += [3*(n.label-1) + ii] if gen_bc[ii]==0 else []
            #Eliminating rows and columns
            K = np.delete( np.delete(K, rows_columns_toNull, 0), rows_columns_toNull, 1 )
            F = np.delete(F ,rows_columns_toNull, 0)
            #Solve Ku=F
            u = np.linalg.inv(K).dot(F)
            #Assembling u again
            uu = np.empty((3 * len(p.nodes)))
            count = 0
            for n in range(3 * len(p.nodes)):
                if any([n == rct for rct in rows_columns_toNull]):
                    uu[n] = 0
                else:
                    uu[n] = u[count]
                    count += 1

            scale_factor = 10
            uu = scale_factor * uu.reshape((len(p.nodes),3))
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            for n in p.nodes:
                ax.scatter(n.coordinates[0] + uu[n.label-1,0], n.coordinates[1] + uu[n.label-1,1],
                           n.coordinates[2] + uu[n.label-1,2], marker='o')



            #Field and history outputs


    return uu


def ElementStiffness(mdlObj, partObj, elementObj):
    u_e = np.array(partObj.nodes[elementObj.connectivity[0]].coordinates)  #Construction of u_e = [x1 y1 z1; x2 y2 z2; x3... ]
    for n in elementObj.connectivity[1:]:
        u_e = np.vstack([u_e , partObj.nodes[n].coordinates])
    # Constant DQ - derivative of the 8 N functions from i, j and k (rows); after replacing i,j,k = 0
    DQ = 0.125 * np.array([[-1, 1, -1, 1, -1, 1, -1, 1],
                           [-1, -1, 1, 1, -1, -1, 1, 1],
                           [-1, -1, -1, -1, 1, 1, 1, 1]])
    #Jacobian inverse already subs zero
    J_inv = np.linalg.inv( DQ.dot(u_e) )
    # Matrics constining the derivatives of N through x,y,z, already with i,j,s = 0
    dQ = np.empty((partObj.elements.nodePerElement,3))
    for n in range(partObj.elements.nodePerElement):
        dQ[n,:] = np.array(np.matrix.transpose(J_inv.dot(DQ[:, n])))

    #Strin-displacemente matrix
    strainDisplacementMatrix = np.array([[dQ[0,0], 0, 0, dQ[1,0], 0, 0, dQ[2,0], 0, 0, dQ[3,0], 0, 0, dQ[4,0], 0, 0, dQ[5,0], 0, 0, dQ[6,0], 0, 0, dQ[7,0], 0, 0],
                      [0, dQ[0,1], 0, 0, dQ[1,1], 0, 0, dQ[2,1], 0, 0, dQ[3,1], 0, 0, dQ[4,1], 0, 0, dQ[5,1], 0, 0, dQ[6,1], 0, 0, dQ[7,1], 0],
                      [0, 0, dQ[0,2], 0, 0, dQ[1,2], 0, 0, dQ[2,2], 0, 0, dQ[3,2], 0, 0, dQ[4,2], 0, 0, dQ[5,2], 0, 0, dQ[6,2], 0, 0, dQ[7,2]],
                      [dQ[0,1], dQ[0,0], 0, dQ[1,1], dQ[1,0], 0, dQ[2,1], dQ[2,0], 0, dQ[3,1], dQ[3,0], 0, dQ[4,1], dQ[4,0], 0, dQ[5,1], dQ[5,0], 0, dQ[6,1], dQ[6,0], 0, dQ[7,1], dQ[7,0], 0],
                      [0, dQ[0,2], dQ[0,1], 0, dQ[1,2], dQ[1,1], 0, dQ[2,2], dQ[2,1], 0, dQ[3,2], dQ[3,1], 0, dQ[4,2], dQ[4,1], 0, dQ[5,2], dQ[5,1], 0, dQ[6,2], dQ[6,1], 0, dQ[7,2], dQ[7,1]],
                      [dQ[0,2], 0, dQ[0,0], dQ[1,2], 0, dQ[1,0], dQ[2,2], 0, dQ[2,0], dQ[3,2], 0, dQ[3,0], dQ[4,2], 0, dQ[4,0], dQ[5,2], 0, dQ[5,0], dQ[6,2], 0, dQ[6,0], dQ[7,2], 0, dQ[7,0]]])

    #Element material
    E = mdlObj.materials[elementObj.material].elastic.table[0]
    v = mdlObj.materials[elementObj.material].elastic.table[1]
    p1 = v * E / ((1 + v) * (1 - 2*v))
    p2 = E / (2 * (1 + v))
    E_e = np.array([[p1 + 2 * p2, p1, p1, 0, 0, 0],
                    [p1, p1 + 2 * p2, p1, 0, 0, 0],
                    [p1, p1, p1 + 2 * p2, 0, 0, 0],
                    [0, 0, 0, p2, 0, 0],
                    [0, 0, 0, 0, p2, 0],
                    [0, 0, 0, 0, 0, p2]])

    K_e = 2*np.matrix.transpose(strainDisplacementMatrix).dot(E_e.dot(strainDisplacementMatrix))
    return K_e



