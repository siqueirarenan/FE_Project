import numpy as np
import matplotlib.pyplot as plt

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
        new_Job = job(job_name)
        return new_Job
    def openMesh(self,file_name):
        self.models[self.models.keys()[0]].acis = None
        pass #Importing step model

class Model:  #class of a model containing all the parts, material properties,... of a mdb
    name = None #name of the model
    parts = {} #dictionary containing all the different parts {'part_name':obj_from_class_part}
    materials = {} #dictionary containing material objects {'material_name':obj_from_class_mat}
    sections = {}
    acis = None #geometry file saved in the model to be imported
    def __init__(self,model_name='Model-1'):
        self.name = model_name
    def Part(self,name='Part-1'):
        self.parts[name] = Part(name)
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

class Part:  #class of parts of a model
    name = None
    geometryFile = None #Input file with geometry of the part
    nodes = None #Object from the class MeshNodeArray
    elements = None #Object from the class MeshElementArray
    sets = {}
    sections ={}
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
        for k in range(depth):
            for j in range(height):
                for i in range(width):
                    self.nodes += [MeshNode((i * ms, j * ms, k * ms), count)]  # Assigning coordinaates
                    count += 1
        self.elements = MeshElementArray([])
        count = 1
        for k in range(depth - 1): #Constructing connectivity vertor and each MeshElement object
            for j in range(height - 1):
                for i in range(width - 1):
                    connectivity = []  #By the iterator, not the label
                    for kk in [0,1]:
                        for jj in [0,1]:
                            for ii in [0,1]:
                                connectivity += [height * width * (k + kk) +
                                                 width * (j + jj) + i + ii]
                    connectivity = tuple(connectivity)
                    self.elements += [MeshElement(connectivity, count)]
                    count += 1
        fig = plt.figure()
        ax = fig.add_subplot('111', projection='3d')
        ax.voxels(np.ones([depth,width,height], dtype=np.bool), facecolors=[1, 1, 1, 1], edgecolors='k')
    def Set(self,set_name,meshElementArrayObj):
        self.sets[set_name] = Set(set_name,meshElementArrayObj)
        return self.sets[set_name]
    def SectionAssignment(self,set_obj,sec_name):
        pass

class Set:
    elements = {}
    def __init__(self,set_name,meshElementArrayObj):
        self.name = set_name
        self.elements = meshElementArrayObj

class MeshElementArray:
    elements = None
    def __init__(self,elements):
        self.elements = elements
    def __getitem__(self, item):
        return self.elements[item]
    def __add__(self, other):
        self.elements += other  #other must be node object!!!
        return self
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
    def getNodes(self): #return a tuble with the node_objects of the element
        pass

class MeshNodeArray:
    nodes = None #Array of MeshNode object
    def __init__(self,nodes):
        self.nodes = nodes
    def __getitem__(self, item):
        return self.nodes[item]
    def __add__(self, other):
        self.nodes += other  #other must be node object!!!
        return self

class MeshNode:
    coordinates = () #Node coordinates as tuple x,y,z
    label = None
    def __init__(self, coordinates, label):
        self.coordinates = coordinates
        self.label = label
    def getElements(self):
        pass

class Section:
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

class job:  #class of jobs of a model
    def __init__(self,job_name):
        self.name = job_name
    def submit(self):  #function that runs the job
        pass

###---------------------------------------------------------------------------------
mdb=Mdb()
# mdb.openMesh("C:/Users/sique/Desktop/PSU_SUP.STL")
# model = mdb['Model-1']
# part = model.PartFromGeometryFile()
# part.seedPart(5)
# part.generateMesh()
mdl=mdb.models['Model-1']
part=mdl.Part('Part-1')
mdl.Material('Material-1').Elastic(((210000,0.3,),))
mdl.HomogeneousSolidSection('Mat2Sec','Material02')
part.uniformHexMesh(10,10,30,1)
part.Set('Set-1',part.elements)




#execfile('C:/Users/sique/Desktop/FE_Project/FE.py')


