import numpy as np
import sympy as sp

#Stiffness matrix construction

def solver(job_name,mdb):

    #Displacement vector of every node in the form [u1x u1y u1z u2x u2y ... ]
    #Initially as the actual coordinate of each node
    #Displacement boundary conditions can be applied

    u = np.array([])
    for p in mdb.models[mdb.jobs[job_name].model].parts.values(): #Construction of u general vector
        for n in p.nodes:
            u = np.append(u,np.array(n.coordinates)) #Linha!

    #Separate the matrix displacement for an element e

    for p in mdb.models[mdb.jobs[job_name].model].parts.values():
        elementssetiffness = {}
        for e in p.elements:                                 #Begining of construction of K_e
            elementssetiffness[e.label] = ElementStiffness(p,e)

    return elementssetiffness


class ElementStiffness:
    def __init__(self,partObj, elementObj):
        self.label = elementObj.label
        # Strain vector s_e = [s_xx s_yy s_zz s_xy s_xz s_yz]
        i, j, k = sp.symbols('i j k')
        Q1 = (1 - i)*(1 - j)*(1 - k) / 8
        Q2 = (1 + i)*(1 - j)*(1 - k) / 8
        Q3 = (1 - i)*(1 + j)*(1 - k) / 8
        Q4 = (1 + i)*(1 + j)*(1 - k) / 8
        Q5 = (1 - i)*(1 - j)*(1 + k) / 8
        Q6 = (1 + i)*(1 - j)*(1 + k) / 8
        Q7 = (1 - i)*(1 + j)*(1 + k) / 8
        Q8 = (1 + i)*(1 + j)*(1 + k) / 8

        N_e = np.array([[Q1, 0, 0, Q2, 0, 0, Q3, 0, 0, Q4, 0, 0, Q5, 0, 0, Q6, 0, 0, Q7, 0, 0, Q8, 0, 0],
                        [0, Q1, 0, 0, Q2, 0, 0, Q3, 0, 0, Q4, 0, 0, Q5, 0, 0, Q6, 0, 0, Q7, 0, 0, Q8, 0],
                        [0, 0, Q1, 0, 0, Q2, 0, 0, Q3, 0, 0, Q4, 0, 0, Q5, 0, 0, Q6, 0, 0, Q7, 0, 0, Q8]])

        self.u_e = np.array([])

        for n in elementObj.connectivity:
            self.u_e = np.append(self.u_e , partObj.nodes[n].coordinates)      #Vector with the coordinates of nodes of the element

        #Constructing the x,y,z in elementwise coordinate system with i,j,k
        self.d_e = sp.Matrix(N_e.dot(self.u_e[:, np.newaxis]))
        #Jacobian
        J = self.d_e.jacobian([i,j,k])
        self.J_inv_easy = J.inv().subs([(i,0),(j,0),(k,0)])
        dQ1 = self.J_inv_easy * sp.Matrix([[sp.diff(Q1, i)], [sp.diff(Q1, j)], [sp.diff(Q1, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ2 = self.J_inv_easy * sp.Matrix([[sp.diff(Q2, i)], [sp.diff(Q1, j)], [sp.diff(Q2, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ3 = self.J_inv_easy * sp.Matrix([[sp.diff(Q3, i)], [sp.diff(Q1, j)], [sp.diff(Q3, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ4 = self.J_inv_easy * sp.Matrix([[sp.diff(Q4, i)], [sp.diff(Q1, j)], [sp.diff(Q4, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ5 = self.J_inv_easy * sp.Matrix([[sp.diff(Q5, i)], [sp.diff(Q1, j)], [sp.diff(Q5, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ6 = self.J_inv_easy * sp.Matrix([[sp.diff(Q6, i)], [sp.diff(Q1, j)], [sp.diff(Q6, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ7 = self.J_inv_easy * sp.Matrix([[sp.diff(Q7, i)], [sp.diff(Q1, j)], [sp.diff(Q7, k)]]).subs([(i,0),(j,0),(k,0)])
        dQ8 = self.J_inv_easy * sp.Matrix([[sp.diff(Q8, i)], [sp.diff(Q1, j)], [sp.diff(Q8, k)]]).subs([(i,0),(j,0),(k,0)])
        #Strin-displacemente matrix
        self.strainDisplacementMatrix = sp.Matrix([[dQ1[0], 0, 0, dQ2[0], 0, 0, dQ3[0], 0, 0, dQ4[0], 0, 0, dQ5[0], 0, 0, dQ6[0], 0, 0, dQ7[0], 0, 0, dQ8[0], 0, 0],
                          [0, dQ1[1], 0, 0, dQ2[1], 0, 0, dQ3[1], 0, 0, dQ4[1], 0, 0, dQ5[1], 0, 0, dQ6[1], 0, 0, dQ7[1], 0, 0, dQ8[1], 0],
                          [0, 0, dQ1[2], 0, 0, dQ2[2], 0, 0, dQ3[2], 0, 0, dQ4[2], 0, 0, dQ5[2], 0, 0, dQ6[2], 0, 0, dQ7[2], 0, 0, dQ8[2]],
                          [dQ1[1], dQ1[0], 0, dQ2[1], dQ2[0], 0, dQ3[1], dQ3[0], 0, dQ4[1], dQ4[0], 0, dQ5[1], dQ5[0], 0, dQ6[1], dQ6[0], 0, dQ7[1], dQ7[0], 0, dQ8[1], dQ8[0], 0],
                          [0, dQ1[2], dQ1[1], 0, dQ2[2], dQ2[1], 0, dQ3[2], dQ3[1], 0, dQ4[2], dQ4[1], 0, dQ5[2], dQ5[1], 0, dQ6[2], dQ6[1], 0, dQ7[2], dQ7[1], 0, dQ8[2], dQ8[1]],
                          [dQ1[2], 0, dQ1[0], dQ2[2], 0, dQ2[0], dQ3[2], 0, dQ3[0], dQ4[2], 0, dQ4[0], dQ5[2], 0, dQ5[0], dQ6[2], 0, dQ6[0], dQ7[2], 0, dQ7[0], dQ8[2], 0, dQ8[0]]])






    #Element shape function
