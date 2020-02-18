import FE_classes, FE_Plots
import time


if __name__=="__main__":
    t_i = time.time()

    mdb = FE_classes.Mdb()

    #Model geometry and properties
    mdl = mdb.models['Model-1']
    part = mdl.Part('Part-1')
    mdl.Material('Material-1').Elastic(table = ((1.0,0.3,),))
    mdl.HomogeneousSolidSection('Section-1','Material-1')
    w,h,d = part.uniformHexMesh(10,10,30,1)
    part.SectionAssignment(part.Set('Set-1',part.elements),'Section-1')

    #Step and outputs
    step = mdl.StaticStep('Step-1')
    mdl.FieldOutputRequest('FieldOut-1', 'Step-1', variables=('MISESMAX','ELEDEN','EVOL', ))

    #Load and boundary conditions
    regL = part.NodeRegionFromFunction(lambda x, y, z: z > -0.001 and z < 0.001 and y > -0.0001 and y < 0.0001)
    mdl.ConcentratedForce('Load-1','Step-1',region = regL,cf1 = 0,cf2 = -1,cf3 = 0)
    regBC = part.NodeRegionFromFunction(lambda x, y, z: z>29.999 and z<30.001)
    mdl.DisplacementBC('BC-1','Step-1',region = regBC,u1=0,u2=0,u3=0) #urs not implemented to make sense, only zero allowed!!!

    #Job
    resuts = mdb.Job('Job-1','Model-1').submit()

    #Plots
    FE_Plots.undeformedNodePlot(mdl, part, step)
    #FE_Plots.deformedNodePlot(mdl, part, step, resuts[0], scale_factor=1)
    FE_Plots.vonMisesHexMeshPlot(w,h,d,resuts[1])


    print(time.time() - t_i)

#execfile('C:/Users/sique/Desktop/FE_Project/FE.py')
