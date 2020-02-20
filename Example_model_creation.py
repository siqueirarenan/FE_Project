import FE_classes, FE_Plots
from FE_classes import openOdb, openMdb
import time

# MODEL CREATION
mdb = FE_classes.Mdb()

# Model geometry and properties
mdl = mdb.models['Model-1']
part = mdl.Part('Part-1')
mdl.Material('Material-1').Elastic(table=((1.0, 0.3,),))
mdl.HomogeneousSolidSection('Section-1', 'Material-1')
part.uniformHexMesh(10, 30, 10, 1)  # MBB beam
part.SectionAssignment(part.Set('Set-1', part.elements), 'Section-1')

# Step and outputs
step = mdl.StaticStep('Step-1')
mdl.FieldOutputRequest('FieldOut-1', 'Step-1', variables=('MISESMAX', 'ESEDEN', 'EVOL', ))
mdl.HistoryOutputRequest('FieldOut-1', 'Step-1', variables=('ALLWK',))

# Load and boundary conditions
regL = part.NodeRegionFromFunction(lambda x, y, z: z == 0.0 and y == 0.0)
mdl.ConcentratedForce('Load-1', 'Step-1', region=regL, cf1=0.0, cf2=0.0, cf3=-1.0)
regBC = part.NodeRegionFromFunction(lambda x, y, z: y == 30.0)
mdl.DisplacementBC('BC-1', 'Step-1', region=regBC, u1=0.0, u2=0.0, u3=0.0) # Only zero allowed

mdb.saveAs('Example_model')
t_i = time.time()

# Job
mdb.Job('Job-1', 'Model-1').submit()
mdb.jobs['Job-1'].waitForCompletion()
odb = openOdb('Output_files\\' + 'Job-1.odb')

# PLOTS
# FE_Plots.undeformedNodePlot(mdl, part, step)
# FE_Plots.deformedNodePlot(mdl, part, step, odb, scale_factor=1)
FE_Plots.FieldOutputHexMeshPlot(part, step, odb, 'MISESMAX')
FE_Plots.FieldOutputHexMeshPlot(part, step, odb, 'ESEDEN')
# FE_Plots.HistoryOutputPlot(odb, step, ['ALLWK'])

print(time.time() - t_i)
