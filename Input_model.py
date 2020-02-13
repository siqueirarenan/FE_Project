import FE_classes
import time


if __name__=="__main__":
    t_i = time.time()

    mdb = FE_classes.Mdb()

    #Model geometry and properties
    mdl = mdb.models['Model-1']
    part = mdl.Part('Part-1')
    mdl.Material('Material-1').Elastic(table = ((210000,0.3,),))
    mdl.HomogeneousSolidSection('Section-1','Material-1')
    part.uniformHexMesh(10,10,30,10)
    part.SectionAssignment(part.Set('Set-1',part.elements),'Section-1')

    #Step and outputs
    mdl.StaticStep('Step-1')
    mdl.FieldOutputRequest('FieldOut-1', 'Step-1', variables=('MISESMAX','ELEDEN','EVOL', ))

    #Load and boundary conditions
    regL = part.NodeRegionFromFunction(lambda x, y, z: z > 29.999 and z < 30.0001 and y > -0.0001 and y < 0.0001)
    mdl.ConcentratedForce('Load-1','Step-1',region = regL,cf1 = 0,cf2 = -1000,cf3 = 0)
    regBC = part.NodeRegionFromFunction(lambda x, y, z: z>-0.0001 and z<0.0001)
    mdl.DisplacementBC('BC-1','Step-1',region = regBC,u1=0,u2=0,u3=0) #urs not implemented to make sense, only zero allowed!!!

    #Job
    u = mdb.Job('Job-1','Model-1').submit()


    print(time.time() - t_i)

#execfile('C:/Users/sique/Desktop/FE_Project/FE.py')