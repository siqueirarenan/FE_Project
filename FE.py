import pymesh
import numpy as np

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
    def openStep(self,file_name):
        self.models[self.models.keys()[0]].acis = pymesh.load_mesh("model.obj") #open(file_name) #Importing step model

class Model:  #class of a model containing all the parts, material properties,... of a mdb
    name = None #name of the model
    parts = {} #dictionary containing all the different parts {'part_name':obj_from_class_part}
    materials = {} #dictionary containing material objects {'material_name':obj_from_class_mat}
    acis = None #geometry file saved in the model to be imported
    def __init__(self,model_name='Model-1'):
        self.name = model_name
    def PartFromGeometryFile(self,name='Part-1'):  #fuction to create a new part
        self.parts[name] = Part(name,self.acis)
        return self.parts[name]
    def Material(self,mat_name): #function to create a new material model
        self.materials[mat_name] = material(mat_name)
        return self.materials[mat_name]

class Part:  #class of parts of a model
    name = None
    geometryFile = None
    nodes = {}
    elements = {}
    sets = {}
    sections ={}
    seed_size = 0
    def __init__(self,name,geometryFile):
        self.name = name
        self.geometry = geometryFile
    def seedPart(self,seed_size):
        self.seed_size = seed_size
    def generateMesh(self):
        box_min = np.array([])
        box_max = box_min + np.array([self.seed_size,self.seed_size,self.seed_size])
        pymesh.generate_box_mesh(box_min, box_max)

class Set(Part):
    pass

class Element:
    pass

class nodes:
    pass

class section:
    pass

class material:  #class of materials of a model
    #como acessar os dados do material
    def __init__(self,mat_name):
        self.name = mat_name
    def Elastic(self,dependencies=1,table=((0, 0,), )):
        pass

class job:  #class of jobs of a model
    def __init__(self,job_name):
        self.name = job_name
    def submit(self):  #function that runs the job
        pass

###---------------------------------------------------------------------------------
mdb=Mdb()
print(mdb.models.items())
model=mdb.models['Model-1']
print(model.name)
part=model.PartFromGeometryFile()
print(part.name)