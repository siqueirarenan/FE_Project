import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #needed
import FE_solver

#mdb = openMdb(path+file) 'C:\Users\sique\Desktop\QLK_Model-614-2.cae'
def openMdb(file_name): #function to open a model data base
    pass

class Mdb:   #class with all the input information in a file (model data base)
    models = {} #dictionary containing the different models in the file {'model_name':obj_from_class_mdl}
    jobs = {}
    def __init__(self):
        self.models['Model-1'] = Model() #iniciation of a mdb creates a standard model
    def saveAs(self,file_name):
        pass
    def Job(self,job_name,model_name): #function to crete a new model
        self.jobs[job_name] = ModelJob(self,job_name,model_name)
        return self.jobs[job_name]
    def openMesh(self,file_name):
        self.models[self.models.keys()[0]].acis = None
        pass #Importing step model

class Model:  #class of a model containing all the parts, material properties,... of a mdb
    name = None #name of the model
    parts = {} #dictionary containing all the different parts {'part_name':obj_from_class_part}
    materials = {} #dictionary containing material objects {'material_name':obj_from_class_mat}
    sections = {}
    steps = {}
    fieldOutputRequests = {}
    historyOutputRequests = {}
    loads ={}
    boundaryConditions = {}
    acis = None #geometry file saved in the model to be imported
    def __init__(self,model_name='Model-1'):
        self.name = model_name
    def Part(self,name='Part-1'):
        self.parts[name] = Part(name)
        self.parts[name].parent = self
        return self.parts[name]
    def PartFromGeometryFile(self,name='Part-1'):  #fuction to create a new part
        self.parts[name] = Part(name,self.acis)
        return self.parts[name]
    def Material(self,mat_name): #function to create a new material model
        self.materials[mat_name] = Material(mat_name)
        return self.materials[mat_name]
    def HomogeneousSolidSection(self,sec_name,mat_name): #function to create a section
        self.sections[sec_name] = Section(sec_name,mat_name)
        return self.sections[sec_name]
    def StaticStep(self,step_name):
        self.steps[step_name] = StaticStep(step_name)
        return self.steps[step_name]
    def FieldOutputRequest(self,fo_name,step_name,variables):
        self.fieldOutputRequests[fo_name] = FieldOutputRequest(fo_name)
        self.steps[step_name]._FieldOutputRequestStates(fo_name,variables)
        return self.fieldOutputRequests[fo_name]
    def HistoryOutputRequest(self,ho_name,step_name,variables):
        self.historyOutputRequests[ho_name] = HistoryOutputRequest(ho_name)
        self.steps[step_name]._HistoryOutputRequestStates(ho_name,variables)
        return self.historyOutputRequests[ho_name]
    def ConcentratedForce(self,name,step_name,region,cf1,cf2,cf3):
        self.loads[name] = ConcentratedForce(name,region)
        self.steps[step_name]._ConcentratedForceState(name,cf1,cf2,cf3)
        return self.loads[name]
    def DisplacementBC(self,name,step_name,region,u1=None,u2=None,u3=None,ur1=None,ur2=None,ur3=None):
        self.boundaryConditions[name] = DisplacementBC(name,region)
        self.steps[step_name]._DisplacementBCState(name,u1,u2,u3,ur1,ur2,ur3)
        return self.boundaryConditions[name]

class Part:  #class of parts of a model
    parent = None
    name = None
    geometryFile = None #Input file with geometry of the part
    nodes = None #Object from the class MeshNodeArray
    elements = None #Object from the class MeshElementArray
    sets = {}
    sectionassignments = []
    seed_size = 0.0 #Seed for the mesh size
    def __init__(self,name,geometryFile=None):
        self.name = name
        self.geometry = geometryFile
    def seedPart(self,seed_size): #Defining mesh size
        self.seed_size = seed_size
    def generateMesh(self):
        mesh = None #pymesh.tetrahedralize(self.geometryFile,self.seed_size) #create a tetraedral mesh
        vert = mesh.vertices
        vox = mesh.voxels
        self.nodes = MeshNodeArray([MeshNode(tuple(vert[i]),i) for i in range(vert.shape[0])])
        self.elements = MeshElementArray([MeshElement(tuple(vox[i]), i) for i in range(vox.shape[0])])
    def uniformHexMesh(self , width , height , depth , ms=1 ):
        # Mesh simple hexahedric uniform
        #ms - mesh size
        width //= ms   #Converting to element quantity unit
        height //= ms
        depth //= ms
        self.nodes = MeshNodeArray([])
        count = 1
        for k in range(depth+1):
            for j in range(height+1):
                for i in range(width+1):
                    self.nodes += [MeshNode((i * ms, j * ms, k * ms), count)]  # Assigning coordinaates
                    count += 1
        self.elements = MeshElementArray([],8)
        count = 1
        for k in range(depth): #Constructing connectivity vertor and each MeshElement object
            for j in range(height):
                for i in range(width):
                    connectivity = []  #By the iterator, not the label
                    for kk in [0,1]:
                        for jj in [0,1]:
                            for ii in [0,1]:
                                connectivity += [(height + 1) * (width + 1) * (k + kk) +
                                                 (width + 1) * (j + jj) + i + ii]
                    connectivity = tuple(connectivity)
                    self.elements += [MeshElement(connectivity, count)]
                    count += 1
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.voxels(np.ones([depth,width,height], dtype=np.bool), facecolors=[1, 1, 1, 1], edgecolors='k')
    def Set(self,set_name,meshElementArrayObj):
        self.sets[set_name] = Set(set_name,meshElementArrayObj,self.name)
        return self.sets[set_name]
    def SetFromElementLabels(self,set_name,label_list):
        self.sets[set_name] = Set(set_name,self.elements.sequenceFromLabels(label_list),self.name)
        return self.sets[set_name]
    def SectionAssignment(self,set_obj,sec_name):
        self.sectionassignments += [SectionAssignment(set_obj,sec_name)]
        for e in set_obj.elements:
            e.setMaterial(self.parent.sections[sec_name].material)
        return self.sectionassignments[-1]
    def NodeRegionFromFunction(self,node_in_function): #create a region of nodes through a boolean function(x,y,z)
        _regionObj = Region(elements=[],nodes=[])
        for i in self.nodes:
            if node_in_function(i.coordinates[0],i.coordinates[1],i.coordinates[2]):
                _regionObj += i
        return _regionObj

class SectionAssignment:
    def __init__(self,set_obj,sec_name):
        self.sectionName = sec_name
        self.region = (set_obj.name,set_obj._parent_part,)

class Set:
    def __init__(self,set_name,meshElementArrayObj,parent_part_name):
        self.name = set_name
        self.elements = meshElementArrayObj
        self._parent_part = parent_part_name

class Region:  #Like a set, but more complex
    elements = None
    nodes = None
    def __init__(self,elements=None,nodes=None):
        self.elements = elements
        self.nodes = nodes
    def __add__(self, other):
        if isinstance(other,MeshNode):
            self.nodes += [other]
        if isinstance(other,MeshElement):
            self.elements += [other]
        return self

class MeshElementArray:
    elements = None
    def __init__(self,elements,npe):
        self.elements = elements
        self.nodePerElement = npe
    def __getitem__(self, item):
        return self.elements[item]
    def __add__(self, other):
        self.elements += other  #other must be element object!!!
        return self
    def __len__(self):
        return len(self.elements)
    def sequenceFromLabels(self,label_list): #retorn a MeshElementArray only with the ones in the label_list
        filtered_ele = list(filter(lambda x:x.label in label_list,self.elements))
        return MeshElementArray(filtered_ele)

class MeshElement:
    connectivity = ()
    label = None
    def __init__(self,connectivity,label):
        self.connectivity = connectivity
        self.label = label
    def getAdjacentElements(self):
        pass
    def setMaterial(self,mat_name):
        self.material = mat_name

class MeshNodeArray:
    nodes = None #Array of MeshNode object
    def __init__(self,nodes):
        self.nodes = nodes
    def __getitem__(self, item):
        return self.nodes[item]
    def __add__(self, other):
        self.nodes += other  #other must be node object!!!
        return self
    def __len__(self):
        return len(self.nodes)

class MeshNode:
    coordinates = () #Node coordinates as tuple x,y,z
    label = None
    def __init__(self, coordinates, label):
        self.coordinates = coordinates
        self.label = label
    def getElements(self):
        pass

class Material:  #class of materials of a model
    name = None
    def __init__(self,mat_name):
        self.name = mat_name
    def Elastic(self,dependencies=1,table=((0, 0,), )):
        self.elastic = Elastic(dependencies,table) #Elastic object
    def Density(self, table=((0, ),)):
        self.density = Density(table)

class Elastic:
    def __init__(self,dependencies,table=((0, 0,), )):
        self.table = table[0] #First value is the E-Modul, second is the poisson

class Density:
    def __init__(self,table=((0,), )):
        self.table = table[0][0]

class Section:
    def __init__(self,sec_name,mat_name):
        self.name = sec_name
        self.material = mat_name

class Step:
    fieldOutputRequestStates = {}
    historyOutputRequestStates = {}
    loadStates = {}
    boundaryConditionStates = {}
    def __init__(self,name):
        self.name = name
    def _FieldOutputRequestStates(self,fo_name,variables):
        self.fieldOutputRequestStates[fo_name] = FieldOutputRequestStates(fo_name,variables)
        return self.fieldOutputRequestStates[fo_name]
    def _HistoryOutputRequestStates(self,ho_name,variables):
        self.historyOutputRequestStates[ho_name] = HistoryOutputRequestStates(ho_name,variables)
        return self.fieldOutputRequestStates[ho_name]
    def _ConcentratedForceState(self,name,cf1,cf2,cf3):
        self.loadStates[name] = ConcentratedForceState(cf1,cf2,cf3)
        return self.loadStates[name]
    def _DisplacementBCState(self,name,u1,u2,u3,ur1,ur2,ur3):
        self.boundaryConditionStates[name] = DisplacementBCState(u1,u2,u3,ur1,ur2,ur3)
        return self.boundaryConditionStates[name]

class StaticStep(Step):
    step_type = 'Static'

class FieldOutputRequest:
    def __init__(self,name):
        self.name = name

class FieldOutputRequestStates:
    def __init__(self,name,variables):
        self.name = name
        self.variables = variables

class HistoryOutputRequest:
    def __init__(self,name):
        self.name = name

class HistoryOutputRequestStates:
    def __init__(self,name,variables):
        self.name = name
        self.variables = variables

class ConcentratedForce:
    def __init__(self,name,region):
        self.name = name
        self.region = region

class ConcentratedForceState:
    def __init__(self,cf1,cf2,cf3):
        self.cf1 = cf1
        self.cf2 = cf2
        self.cf3 = cf3

class DisplacementBC:
    def __init__(self,name,region):
        self.name = name
        self.region = region

class DisplacementBCState:
    def __init__(self,u1=None,u2=None,u3=None,ur1=None,ur2=None,ur3=None): #Empty values mean no BC set
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3
        self.ur1 = ur1
        self.ur2 = ur2
        self.ur3 = ur3

class ModelJob:  #class of jobs of a model
    def __init__(self,mdb_obj,job_name,model_name):
        self.name = job_name
        self.model = model_name
        self.mdb = mdb_obj
    def submit(self):  #function that runs the job
        elementsstiffness,K = FE_solver.solver(self.name, self.mdb)
        return u


