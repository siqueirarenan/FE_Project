import numpy as np
import FE_classes
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import sympy as sp
from scipy import linalg as spln
import Gauss

#Stiffness matrix construction

def solver(job_name,mdb):
    mdl = mdb.models[mdb.jobs[job_name].model]
    for p in mdl.parts.values(): #Construction of u general vector
        #Stiffness Matrix construction
        elementsstiffness = {}
        for e in p.elements:                                 #Begining of construction of K_e
            elementsstiffness[e.label] = ElementStiffness(mdl,p,e)
        #Assembly
        K = np.zeros((3*len(p.nodes),3*len(p.nodes)))
        for k_e in elementsstiffness.items():
            for m in range(p.elements.nodePerElement):
                for n in range(p.elements.nodePerElement): #Only up half of the simmetrtic matrix is calculated
                    K[3 * p.elements[k_e[0]-1].connectivity[m]:3 * p.elements[k_e[0]-1].connectivity[m] + 3,
                      3 * p.elements[k_e[0]-1].connectivity[n]:3 * p.elements[k_e[0]-1].connectivity[n] + 3
                    ] += k_e[1][3 * m:3 * m + 3, 3 * n :3 * n + 3]
        #Step
        for s in mdl.steps.values(): #per Step
            #Load vector
            F = np.zeros((3 * len(p.nodes),1))
            for f in s.loadStates.items(): #Checking each load in the step
                if isinstance(f[1],FE_classes.ConcentratedForceState): #Dealing according to case
                    gen_force = np.array([f[1].cf1,f[1].cf2 ,f[1].cf3]) / len(
                        mdl.loads[f[0]].region.nodes)
                for n in mdl.loads[f[0]].region.nodes:  #checking each node of the load
                    F[3 * (n.label-1):3 * (n.label-1) + 3 ,0] +=  gen_force

            #Boundary conditions
            #u = np.empty((3 * len(p.nodes), 1))
            rows_columns_toNull = []
            for bc in s.boundaryConditionStates.items():  # Checking each BC in the step
                if isinstance(bc[1], FE_classes.DisplacementBCState):  # Dealing according to case
                    gen_bc = np.array([bc[1].u1,bc[1].u2 ,bc[1].u3])
                for n in mdl.boundaryConditions[bc[0]].region.nodes:  #checking each node of the load
                    #u[3 * (n.label-1):3 * (n.label-1) + 3 ,0] +=  gen_bc #still not implemented for values != 0
                    #getting the doF with null displacement
                    for ii in range(3):
                        rows_columns_toNull += [3*(n.label-1) + ii] if gen_bc[ii]==0 else []
            #Eliminating rows and columns
            K = np.delete( np.delete(K, rows_columns_toNull, 0), rows_columns_toNull, 1 )
            F = np.delete(F ,rows_columns_toNull, 0)
            #for r in rows_columns_toNull:
            #    K[r,r] = 1e36

            #SOLVER Ku=F
            u = spln.solve(K , F , assume_a='sym')

            #Assembling u again
            uu = np.zeros((3 * len(p.nodes)))
            count = 0
            for dof in range(3 * len(p.nodes)):
                if any([dof == rct for rct in rows_columns_toNull]):
                    uu[dof] = 0
                else:
                    uu[dof] = u[count]
                    count += 1
            #uu=u

            #PLOTS (do in other file!!!)
            #Undeformed nodes showing loads and BC
            fig = plt.figure()
            ax = fig.add_subplot(121, projection='3d')
            for n in p.nodes: #Ploting nodes
                ax.scatter(n.coordinates[0], n.coordinates[1],
                           n.coordinates[2], marker='.', color=[0.3,0.3,0.3])
            for l in s.loadStates.items():  #Marking load nodes
                for n in mdl.loads[l[0]].region.nodes:
                    ax.scatter(n.coordinates[0], n.coordinates[1],
                               n.coordinates[2], marker='o', color='blue')
            for bc in s.boundaryConditionStates.items():  #Marking BC nodes
                for n in mdl.boundaryConditions[bc[0]].region.nodes:
                    ax.scatter(n.coordinates[0], n.coordinates[1],
                               n.coordinates[2], marker='o', color='red')
            #Deformed nodes
            scale_factor = 1
            uu = uu.reshape((len(p.nodes),3))
            ug = scale_factor * uu
            ax = fig.add_subplot(122, projection='3d')
            n_size_o = []
            for n in p.nodes: #Calculate the norm of the displacements
                n_size_o += [np.linalg.norm(uu[n.label - 1,:])]
            n_size = n_size_o / max(n_size_o)
            for n in p.nodes:
                ax.scatter(n.coordinates[0] + ug[n.label-1,0], n.coordinates[1] + ug[n.label-1,1],
                           n.coordinates[2] + ug[n.label-1,2], marker='.', color=np.array([n_size[n.label-1],0,1-n_size[n.label-1]]))

            for l in s.loadStates.items():  #Marking load nodes
                for n in mdl.loads[l[0]].region.nodes:
                    ax.scatter(n.coordinates[0] + ug[n.label-1,0], n.coordinates[1] + ug[n.label-1,1],
                               n.coordinates[2] + ug[n.label-1,2], marker='o', color='blue')
            for bc in s.boundaryConditionStates.items():  #Marking BC nodes
                for n in mdl.boundaryConditions[bc[0]].region.nodes:
                    ax.scatter(n.coordinates[0] + ug[n.label-1,0], n.coordinates[1] + ug[n.label-1,1],
                               n.coordinates[2] + ug[n.label-1,2], marker='o', color='red')


            #Field and history outputs


    return (elementsstiffness,K,F,u,uu,n_size_o)


def myInverse(A):
    detA = A[0, 0] * A[1, 1] * A[2, 2] + A[1, 0] * A[2, 1] * A[0, 2] + A[2, 0] * A[0, 1] * A[1, 2] - A[0,
                0] * A[2, 1] * A[1, 2] - A[2, 0] * A[1, 1] * A[0, 2] - A[1, 0] * A[0, 1] * A[2, 2]
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

from numpy import linalg as la

def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3

def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False

def ElementStiffness(mdlObj, partObj, elementObj):
    i,j,k = sp.symbols('i j k')
    N={}
    N[0] = (1-i)*(1-j)*(1-k)/8
    N[1] = (1+i)*(1-j)*(1-k)/8
    N[2] = (1-i)*(1+j)*(1-k)/8
    N[3] = (1+i)*(1+j)*(1-k)/8
    N[4] = (1-i)*(1-j)*(1+k)/8
    N[5] = (1+i)*(1-j)*(1+k)/8
    N[6] = (1-i)*(1+j)*(1+k)/8
    N[7] = (1+i)*(1+j)*(1+k)/8
    NN = N[0]*np.eye(3)
    for n in range(1,8):
        NN = np.concatenate((NN,N[n]*np.eye(3)),axis=1)
    #Nodes coordenates
    c_xyz = np.array(partObj.nodes[elementObj.connectivity[0]].coordinates)  # Construction of u_e = [x1 y1 z1; x2 y2 z2; x3... ]
    for n in elementObj.connectivity[1:]:
        c_xyz = np.hstack([c_xyz, partObj.nodes[n].coordinates])
    c_ijk = NN.dot(c_xyz)
    #Jacobian
    J = np.array([[sp.diff(c_ijk[0],i) , sp.diff(c_ijk[1],i) , sp.diff(c_ijk[2],i)],
                  [sp.diff(c_ijk[0],j) , sp.diff(c_ijk[1],j) , sp.diff(c_ijk[2],j)],
                  [sp.diff(c_ijk[0],k) , sp.diff(c_ijk[1],k) , sp.diff(c_ijk[2],k)]])
    J_inv, detJ = myInverse(J)
    dQ=[]
    #Derivatives
    for n in range(8):
        dNn = np.array([[sp.diff(N[n],i)],[sp.diff(N[n],j)],[sp.diff(N[n],k)]])
        dNn_xyz = np.dot(J_inv,dNn)
        dQ += [dNn_xyz]
    B = np.array([[dQ[0][0], 0, 0, dQ[1][0], 0, 0, dQ[2][0], 0, 0, dQ[3][0], 0, 0, dQ[4][0], 0, 0, dQ[5][0], 0, 0, dQ[6][0], 0, 0, dQ[7][0], 0, 0],
                  [0, dQ[0][1], 0, 0, dQ[1][1], 0, 0, dQ[2][1], 0, 0, dQ[3][1], 0, 0, dQ[4][1], 0, 0, dQ[5][1], 0, 0, dQ[6][1], 0, 0, dQ[7][1], 0],
                  [0, 0, dQ[0][2], 0, 0, dQ[1][2], 0, 0, dQ[2][2], 0, 0, dQ[3][2], 0, 0, dQ[4][2], 0, 0, dQ[5][2], 0, 0, dQ[6][2], 0, 0, dQ[7][2]],
                  [0, dQ[0][2], dQ[0][1], 0, dQ[1][2], dQ[1][1], 0, dQ[2][2], dQ[2][1], 0, dQ[3][2], dQ[3][1], 0, dQ[4][2], dQ[4][1], 0, dQ[5][2], dQ[5][1], 0, dQ[6][2], dQ[6][1], 0, dQ[7][2], dQ[7][1]],
                  [dQ[0][2], 0, dQ[0][0], dQ[1][2], 0, dQ[1][0], dQ[2][2], 0, dQ[2][0], dQ[3][2], 0, dQ[3][0], dQ[4][2], 0, dQ[4][0], dQ[5][2], 0, dQ[5][0], dQ[6][2], 0, dQ[6][0], dQ[7][2], 0, dQ[7][0]],
                  [dQ[0][1], dQ[0][0], 0, dQ[1][1], dQ[1][0], 0, dQ[2][1], dQ[2][0], 0, dQ[3][1],dQ[3][0], 0, dQ[4][1], dQ[4][0], 0, dQ[5][1], dQ[5][0], 0, dQ[6][1], dQ[6][0], 0, dQ[7][1], dQ[7][0], 0]])
    #Integral
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

    K_ni = np.dot(np.dot(np.matrix.transpose(B),np.dot(E_e,B)),detJ)
    K = np.zeros((K_ni.shape[0],K_ni.shape[1]))
    for n1 in range(K_ni.shape[0]):
        for n2 in range(K_ni.shape[1]):
            ## 3 Gauss points (25s for 1 element)
            #p1 = (8/9)*K_ni[n1,n2][0].subs(i,0) + (5/9)*K_ni[n1,n2][0].subs(i,(3/5)**(0.5)) + (
            #        5/9)*K_ni[n1,n2][0].subs(i,-(3/5)**(0.5))
            #p2 = (8/9)*p1.subs(j,0) + (5/9)*p1.subs(j,(3/5)**(0.5)) + (5/9)*p1.subs(j,-(3/5)**(0.5))
            #K[n1,n2] += [(8/9)*p2.subs(k,0) + (5/9)*p2.subs(k,(3/5)**(0.5)) + (5/9)*p2.subs(k,-(3/5)**(0.5))]
            ## 2 Gauss points (15s)
            #p1 = (K_ni[n1,n2][0].subs(i,-1/(3**0.5)) + K_ni[n1,n2][0].subs(i,1/(3**0.5)))
            #p2 = (p1.subs(j,-1/(3**0.5)) + p1.subs(j,1/(3**0.5)))
            #K[n1,n2] += [(p2.subs(k, -1/(3**0.5)) + p2.subs(k, 1/(3**0.5)))]
            ## 1 Gauss point (3s)
            p1 = 2*K_ni[n1,n2][0].subs(i,0)
            p2 = 2*p1.subs(j,0)
            K[n1,n2] += [2*p2.subs(k,0)]

            K = nearestPD(K)

    return K
