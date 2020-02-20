import FE_solver
import pickle
from tkinter import Tk, Label, Entry, Button, W


def openMdb(file_name):  # function to open a model data base
    with open(file_name, 'rb') as mdbinput:
        mdb = pickle.load(mdbinput)
    return mdb


def openOdb(file_name):
    with open(file_name, 'rb') as odbinput:
        odb = pickle.load(odbinput)
    return odb


def getInputs(parameters, dialogTitle='Inputs'):
    window = Tk()
    window.title('BESO for FE Project')
    r = 0
    entries = []
    for p in parameters:
        Label(window, text=p[0]).grid(row=r, column=0)
        e = Entry(window)
        e.insert(10, p[1])
        e.grid(row=r, column=1, sticky=W)
        entries += [e.get()]
        r += 1
    OKbutton = Button(window, text='OK', command=window.quit)
    OKbutton.grid(row=r, column=0, sticky=W)
    window.mainloop()
    window.destroy()
    return entries


class Mdb:  # class with all the input information in a file (model data base)
    def __init__(self):
        self.models = {'Model-1': Model()}  # iniciation of a mdb creates a standard model
        self.jobs = {}

    def saveAs(self, file_name):
        with open(file_name + '.rcae', 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def Job(self, job_name, model_name):  # function to crete a new model
        self.jobs[job_name] = ModelJob(self, job_name, model_name)
        return self.jobs[job_name]


class Model:  # class of a model containing all the parts, material properties,... of a mdb
    def __init__(self, model_name='Model-1'):
        self.name = model_name
        self.parts = {}  # dictionary containing all the different parts {'part_name':obj_from_class_part}
        self.materials = {}  # dictionary containing material objects {'material_name':obj_from_class_mat}
        self.sections = {}
        self.steps = {}
        self.fieldOutputRequests = {}
        self.historyOutputRequests = {}
        self.loads = {}
        self.boundaryConditions = {}

    def Part(self, name='Part-1'):
        self.parts[name] = Part(name)
        self.parts[name].parent = self
        return self.parts[name]

    def Material(self, mat_name):  # function to create a new material model
        self.materials[mat_name] = Material(mat_name)
        return self.materials[mat_name]

    def HomogeneousSolidSection(self, sec_name, mat_name):  # function to create a section
        self.sections[sec_name] = Section(sec_name, mat_name)
        return self.sections[sec_name]

    def StaticStep(self, step_name):
        self.steps[step_name] = StaticStep(step_name)
        return self.steps[step_name]

    def FieldOutputRequest(self, fo_name, step_name, variables):
        self.fieldOutputRequests[fo_name] = FieldOutputRequest(fo_name)
        self.steps[step_name]._FieldOutputRequestStates(fo_name, variables)
        return self.fieldOutputRequests[fo_name]

    def HistoryOutputRequest(self, ho_name, step_name, variables):
        self.historyOutputRequests[ho_name] = HistoryOutputRequest(ho_name)
        self.steps[step_name]._HistoryOutputRequestStates(ho_name, variables)
        return self.historyOutputRequests[ho_name]

    def ConcentratedForce(self, name, step_name, region, cf1, cf2, cf3):
        self.loads[name] = ConcentratedForce(name, region)
        self.steps[step_name]._ConcentratedForceState(name, cf1, cf2, cf3)
        return self.loads[name]

    def DisplacementBC(self, name, step_name, region, u1=None, u2=None, u3=None, ur1=None, ur2=None, ur3=None):
        self.boundaryConditions[name] = DisplacementBC(name, region)
        self.steps[step_name]._DisplacementBCState(name, u1, u2, u3, ur1, ur2, ur3)
        return self.boundaryConditions[name]


class Part:  # class of parts of a model
    def __init__(self, name, geometryFile=None):
        self.name = name
        self.geometry = geometryFile
        self.parent = None
        self.geometryFile = None  # Input file with geometry of the part
        self.nodes = None  # Object from the class MeshNodeArray
        self.elements = None  # Object from the class MeshElementArray
        self.sets = {}
        self.sectionassignments = []
        self.elementType = None
        self.seed_size = None
        self.HexElementVolume = None
        self.HexElementSizes = None

    def seedPart(self, seed_size):  # Defining mesh size
        self.seed_size = seed_size

    def generateMesh(self):
        mesh = None  # pymesh.tetrahedralize(self.geometryFile,self.seed_size) #create a tetraedral mesh
        vert = mesh.vertices
        vox = mesh.voxels
        self.nodes = MeshNodeArray([MeshNode(tuple(vert[i]), i) for i in range(vert.shape[0])])
        self.elements = MeshElementArray([MeshElement(tuple(vox[i]), i) for i in range(vox.shape[0])])

    def uniformHexMesh(self, w, h, d, ms):
        self.elementType = 'uniformHex'
        # Mesh simple hexahedric uniform
        width = round(w / ms)  # Calculating number of elements on each side
        height = round(h / ms)
        depth = round(d / ms)
        ms_w, ms_h, ms_d = w / width, h / height, d / depth  # calculating element sizes
        self.HexElementVolume = ms_w * ms_h * ms_d
        self.HexElementSizes = (width, height, depth)
        # print(ms_w, ms_h, ms_d)
        self.nodes = MeshNodeArray([])
        count = 1
        for k in range(depth + 1):
            for j in range(height + 1):
                for i in range(width + 1):
                    self.nodes += [MeshNode((i * ms_w, j * ms_h, k * ms_d), count)]  # Assigning coordinates
                    count += 1
        self.elements = MeshElementArray([], 8)
        count = 1
        for k in range(depth):  # Constructing connectivity vertor and each MeshElement object
            for j in range(height):
                for i in range(width):
                    connectivity = []  # By the iterator, not the label
                    for kk in [0, 1]:
                        for jj in [0, 1]:
                            for ii in [0, 1]:
                                connectivity += [(height + 1) * (width + 1) * (k + kk) +
                                                 (width + 1) * (j + jj) + i + ii]
                    connectivity = tuple(connectivity)
                    self.elements += [MeshElement(connectivity, count)]
                    count += 1
        print('{} hexaedric elements created'.format(count - 1))
        return width, height, depth

    def Set(self, set_name, meshElementArrayObj):
        self.sets[set_name] = Set(set_name, meshElementArrayObj, self.name)
        return self.sets[set_name]

    def SetFromElementLabels(self, set_name, label_list):
        self.sets[set_name] = Set(set_name, self.elements.sequenceFromLabels(label_list), self.name)
        return self.sets[set_name]

    def SectionAssignment(self, set_obj, sec_name):
        self.sectionassignments += [SectionAssignment(set_obj, sec_name)]
        for e in set_obj.elements:
            e.setMaterial(self.parent.sections[sec_name].material)
        return self.sectionassignments[-1]

    def NodeRegionFromFunction(self, node_in_function):  # create a region of nodes through a boolean function(x,y,z)
        _regionObj = Region(elements=[], nodes=[])
        for i in self.nodes:
            if node_in_function(i.coordinates[0], i.coordinates[1], i.coordinates[2]):
                _regionObj += i
        return _regionObj


class SectionAssignment:
    def __init__(self, set_obj, sec_name):
        self.sectionName = sec_name
        self.region = (set_obj.name, set_obj._parent_part,)


class Set:
    def __init__(self, set_name, meshElementArrayObj, parent_part_name):
        self.name = set_name
        self.elements = meshElementArrayObj
        self._parent_part = parent_part_name


class Region:  # Like a set, but more complex
    def __init__(self, elements=None, nodes=None):
        self.elements = elements
        self.nodes = nodes

    def __add__(self, other):
        if isinstance(other, MeshNode):
            self.nodes += [other]
        if isinstance(other, MeshElement):
            self.elements += [other]
        return self


class MeshElementArray:
    def __init__(self, elements=None, npe=None):
        self.elements = elements
        self.nodePerElement = npe

    def __getitem__(self, item):
        return self.elements[item]

    def __add__(self, other):
        assert (isinstance(other[0], MeshElement))
        self.elements += other
        return self

    def __len__(self):
        return len(self.elements)

    def sequenceFromLabels(self, label_list):  # retorn a MeshElementArray only with the ones in the label_list
        filtered_ele = list(filter(lambda x: x.label in label_list, self.elements))
        return MeshElementArray(filtered_ele)


class MeshElement:
    _objs = set()

    def __init__(self, connectivity, label):
        self.connectivity = connectivity
        self.label = label
        self._objs.add(self)
        self.material = None

    def setMaterial(self, mat_name):
        self.material = mat_name

    def getAdjacentElements(self):
        adj_e = []
        for e in self._objs:
            if len(set(self.connectivity).intersection(set(e.connectivity))) >= 3 and e is not self:
                adj_e += [e]
        return MeshElementArray(adj_e, npe=8)


class MeshNodeArray:
    def __init__(self, nodes):
        self.nodes = nodes

    def __getitem__(self, item):
        return self.nodes[item]

    def __add__(self, other):
        assert (isinstance(other[0], MeshNode))
        self.nodes += other
        return self

    def __len__(self):
        return len(self.nodes)


class MeshNode:
    def __init__(self, coordinates, label):
        self.coordinates = coordinates
        self.label = label

    def getElements(self):
        pass


class Material:  # class of materials of a model
    def __init__(self, mat_name):
        self.name = mat_name
        self.elastic = None
        self.density = None

    def Elastic(self, table=((0, 0,),), dependencies=1):
        self.elastic = Elastic(dependencies, table)  # Elastic object

    def Density(self, table=((0,),)):
        self.density = Density(table)


class Elastic:
    def __init__(self, dependencies=None, table=((0, 0,),)):
        self.table = table[0]  # First value is the E-Modul, second is the poisson
        self.dependencies = dependencies


class Density:
    def __init__(self, table=((0,),)):
        self.table = table[0][0]


class Section:
    def __init__(self, sec_name, mat_name):
        self.name = sec_name
        self.material = mat_name


class Step:
    def __init__(self, name):
        self.name = name
        self.fieldOutputRequestStates = {}
        self.historyOutputRequestStates = {}
        self.loadStates = {}
        self.boundaryConditionStates = {}

    def _FieldOutputRequestStates(self, fo_name, variables):
        self.fieldOutputRequestStates[fo_name] = FieldOutputRequestStates(fo_name, variables)
        return self.fieldOutputRequestStates[fo_name]

    def _HistoryOutputRequestStates(self, ho_name, variables):
        self.historyOutputRequestStates[ho_name] = HistoryOutputRequestStates(ho_name, variables)
        return self.historyOutputRequestStates[ho_name]

    def _ConcentratedForceState(self, name, cf1, cf2, cf3):
        self.loadStates[name] = ConcentratedForceState(cf1, cf2, cf3)
        return self.loadStates[name]

    def _DisplacementBCState(self, name, u1, u2, u3, ur1, ur2, ur3):
        self.boundaryConditionStates[name] = DisplacementBCState(u1, u2, u3, ur1, ur2, ur3)
        return self.boundaryConditionStates[name]


class StaticStep(Step):
    step_type = 'Static'


class FieldOutputRequest:
    def __init__(self, name):
        self.name = name


class FieldOutputRequestStates:
    def __init__(self, name, variables):
        self.name = name
        self.variables = variables


class HistoryOutputRequest:
    def __init__(self, name):
        self.name = name


class HistoryOutputRequestStates:
    def __init__(self, name, variables):
        self.name = name
        self.variables = variables


class ConcentratedForce:
    def __init__(self, name, region):
        self.name = name
        self.region = region


class ConcentratedForceState:
    def __init__(self, cf1, cf2, cf3):
        self.cf1 = cf1
        self.cf2 = cf2
        self.cf3 = cf3


class DisplacementBC:
    def __init__(self, name, region):
        self.name = name
        self.region = region


class DisplacementBCState:
    def __init__(self, u1=None, u2=None, u3=None, ur1=None, ur2=None, ur3=None):  # Empty values mean no BC set
        self.u1 = u1
        self.u2 = u2
        self.u3 = u3
        self.ur1 = ur1
        self.ur2 = ur2
        self.ur3 = ur3


class ModelJob:  # class of jobs of a model
    def __init__(self, mdb_obj, job_name, model_name):
        self.name = job_name
        self.model = model_name
        self.mdb = mdb_obj

    def submit(self):  # function that runs the job
        FE_solver.solver(self.name, self.mdb)

    def waitForCompletion(self):
        pass


class Odb:
    def __init__(self, name):
        self.steps = {}
        self.name = name
        self.rootAssembly = OdbAssembly()

    def Step(self, name):
        self.steps[name] = OdbStep(name)
        return self.steps[name]

    def save(self):
        with open(self.name + '.odb', 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def close(self):
        f = open(self.name + '.odb', 'rb')
        f.close()


class OdbStep:
    def __init__(self, name):
        self.frames = []
        self.historyRegions = {}
        self.name = name

    def Frame(self):
        frame_number = len(self.frames)
        self.frames += [OdbFrame(frame_number)]
        return self.frames[-1]

    def HistoryRegion(self, name):
        self.historyRegions[name] = HistoryRegion(name)
        return self.historyRegions[name]


SCALAR = None


class OdbFrame:
    def __init__(self, frame_number):
        self.number = frame_number
        self.fieldOutputs = {}

    def FieldOutput(self, name, description=None, fotype=SCALAR):
        self.fieldOutputs[name] = FieldOutput(name, description, fotype)
        return self.fieldOutputs[name]


WHOLE_ELEMENT = None


class FieldOutput:
    def __init__(self, name, description, fotype):
        self.name = name
        self.values = FieldValueArray()
        self.description = description
        self.fotype = fotype

    def addData(self, position, instance, labels, data):
        self.values.ElementArray(data, labels)


class FieldValueArray:
    def __init__(self):
        self.original_data = None
        self.fieldvalues = []
        self.fieldvalues = []

    def NodeArray(self, orderedfieldvalues, labels=0):
        self.original_data = orderedfieldvalues
        if labels == 0:
            labels = list(range(1, 1 + len(orderedfieldvalues)))
        count = 0
        for i in orderedfieldvalues:
            self.fieldvalues += [FieldValue(labels[count]).NodeData(i)]
            count += 1
        return self

    def ElementArray(self, orderedfieldvalues, labels=0):
        self.original_data = orderedfieldvalues
        if labels == 0:
            labels = list(range(1, 1 + len(orderedfieldvalues)))
        count = 0
        for i in orderedfieldvalues:
            self.fieldvalues += [FieldValue(labels[count]).ElementData(i)]
            count += 1
        return self

    def __getitem__(self, item):
        return self.fieldvalues[item]


class FieldValue:
    def __init__(self, name):
        self.name = name
        self.data = None
        self.nodeLabel = None
        self.elementLabel = None

    def NodeData(self, data):
        self.data = data
        self.nodeLabel = self.name
        return self

    def ElementData(self, data):
        self.data = data
        self.elementLabel = self.name
        return self


class OdbAssembly:
    def __init__(self):
        self.instances = {}
        self.elements = {}
        self.nodes = {}

    def Instance(self, name, obj=None):
        self.instances[name] = OdbInstance(name, obj)
        return self.instances[name]


class OdbInstance:
    def __init__(self, name, obj):
        self.name = name
        self.obj = obj
        self.elementSets = {}
        self.nodeSets = {}

    def ElementSet(self, name, meshArrayObj):
        self.elementSets[name] = OdbSet(name, meshArrayObj)

    def NodeSet(self, name, meshArrayObj):
        self.nodeSets[name] = OdbSet(name, meshArrayObj)


class OdbSet:
    def __init__(self, name, meshArrayObj):
        self.name = name
        self.elements = meshArrayObj if isinstance(meshArrayObj, MeshElementArray) else None
        self.nodes = meshArrayObj if isinstance(meshArrayObj, MeshNodeArray) else None


class HistoryRegion:
    def __init__(self, name):
        self.name = name
        self.historyOutputs = {}

    def HistoryOutput(self, name, frame, data):
        self.historyOutputs[name] = HistoryOutput(name, frame, data)
        return self.historyOutputs[name]


class HistoryOutput:
    def __init__(self, name, frame=0, data=None):
        self.name = name
        self.data = ()
        self.addData(frame, data)

    def addData(self, frame, data):
        data_list = list(self.data)
        data_list += [(frame, data,)]
        self.data = tuple(data_list)
