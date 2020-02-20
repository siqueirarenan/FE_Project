"""General 3D topology optimization code by Zhihao Zuo and Yimin Xie. Note that the CAE file shall contain a model 'Model-1' with a dependent part 'Part-1' and a static step 'Step-1'.
Lightly adapted version for FE Project, by Renan Siqueira"""
import math, sys
from FE_classes import *
from FE_Plots import *
import matplotlib.pyplot as plt
## Function of formatting Abaqus model for stiffness optimisation
def fmtMdb(Mdb):
    mdl = Mdb.models['Model-1']
    part = mdl.parts['Part-1']
    # Build sections and assign solid section
    mdl.Material('Material01').Elastic(((1.0, 0.3), ))
    mdl.HomogeneousSolidSection('sldSec','Material01')
    mdl.Material('Material02').Elastic(((0.001**3, 0.3), ))
    mdl.HomogeneousSolidSection('voidSec','Material02')
    part.SectionAssignment(part.Set('ss',part.elements),'sldSec')
    # Define output request
    mdl.FieldOutputRequest('SEDensity','Step-1',variables=('ESEDEN', ))
    mdl.HistoryOutputRequest('ExtWork','Step-1',variables=('ALLWK', ))
## Function of running FEA for raw sensitivities and objective function
def FEA(Iter,Mdb,Xe,Ae):
    Mdb.Job('Design_Job'+str(Iter),'Model-1').submit()
    Mdb.jobs['Design_Job'+str(Iter)].waitForCompletion()
    opdb = openOdb('Output_files\\' + 'Design_Job'+str(Iter)+'.odb')
    seng = opdb.steps['STEP-1'].frames[-1].fieldOutputs['ESEDEN'].values
    for en in seng: Ae[en.elementLabel]=en.data/Xe[en.elementLabel]
    obj=opdb.steps['STEP-1'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLWK'].data[-1][1]
    opdb.close()
    return obj
## Function of preparing filter map (Fm={elm1:[[el1,el2,...],[wf1,wf2,...]],...})
def preFlt(Rmin,Elmts,Nds,Fm):
    # Calculate element centre coordinates
    c0 = {}
    for el in Elmts:
        nds = el.connectivity
        c0[el.label]=[sum([Nds[nd].coordinates[i]/len(nds) for nd in nds]) for i in range(3)]
    # Weighting factors
    for el in Elmts:
        Fm[el.label] = [[],[]]
        for em in Elmts:
            dis=math.sqrt(sum([(c0[el.label][i]-c0[em.label][i])**2 for i in range(3)]))
            if dis<Rmin:
                Fm[el.label][0].append(em.label)
                Fm[el.label][1].append(Rmin - dis)
        sm = sum(Fm[el.label][1])
        for i in range(len(Fm[el.label][0])): Fm[el.label][1][i] /= sm
## Function of filtering sensitivities
def fltAe(Ae,Fm):
    raw = Ae.copy()
    for el in Fm.keys():
        Ae[el] = 0.0
        for i in range(len(Fm[el][0])): Ae[el]+=raw[Fm[el][0][i]]*Fm[el][1][i]
## Function of optimality update for design variables and Abaqus model
def BESO(Vf,Xe,Ae,Part,Elmts):
    lo, hi = min(Ae.values()), max(Ae.values())
    tv = Vf*len(Elmts)
    while (hi-lo)/hi > 1.0e-5:
        th = (lo+hi)/2.0
        for key in Xe.keys(): Xe[key] = 1.0 if Ae[key]>th else 0.001
        if sum(Xe.values())-tv>0: lo = th
        else: hi = th
    # Label elements as solid or void
    vlb, slb = [], []
    for el in Elmts:
        if Xe[el.label] == 1.0: slb.append(el.label)
        else: vlb.append(el.label)
    # Assign solid and void elements to each section
    Part.SectionAssignment(Part.SetFromElementLabels('ss',slb),'sldSec')
    Part.SectionAssignment(Part.SetFromElementLabels('vs',vlb),'voidSec')
## ====== MAIN PROGRAM ======
if __name__ == '__main__':
    # Set parameters and inputs
    pars = (('VolFrac:','0.5'), ('Rmin:', '1'), ('ER:', '0.02'))
    vf,rmin,ert = [float(k) if k!=None else 0 for k in getInputs(pars,dialogTitle='Parameters')]
    if vf<=0 or rmin<0 or ert<=0: sys.exit("Bad input")
    mddb = openMdb(getInputs((('Input CAE file:','C:\\Users\\sique\\Desktop\\FE_Project\\Example_model.rcae'),))[0])
    # Design initialization
    fmtMdb(mddb)
    part = mddb.models['Model-1'].parts['Part-1']
    elmts, nds = part.elements, part.nodes
    oh, vh = [], []
    xe, ae, oae, fm = {}, {}, {}, {}
    for el in elmts: xe[el.label] = 1.0
    if rmin>0: preFlt(rmin,elmts,nds,fm)
    # Optimisation iteration
    change, iter, obj = 1, -1, 0
    geometry_fig = plt.figure(figsize=[10,10])
    while change > 0.001:
        iter += 1
        # Run FEA
        oh.append(FEA(iter,mddb,xe,ae))
        # Process sensitivities
        if rmin>0: fltAe(ae,fm)
        if iter > 0: ae=dict([(k,(ae[k]+oae[k])/2.0) for k in ae.keys()])
        oae = ae.copy()
        # BESO optimisation
        vh.append(sum(xe.values())/len(xe))
        nv = max(vf,vh[-1]*(1.0-ert))
        BESO(nv,xe,ae,part,elmts)
        if iter>10: change=math.fabs((sum(oh[iter-4:iter+1])-sum(oh[iter-9:iter-4]))/sum(oh[iter-9:iter-4]))
        print('Iteration: {}  -  Volume: {}  -  Compliance: {}'.format(iter,vh[-1],oh[-1]))
        HexGeometryPlot(geometry_fig, part, ['ss'])
        plt.savefig('Output_files\\' + 'Iteration_' + str(iter) + '.png',bbox_inches = 'tight',)
    # Save results
    with open('Output_files\\' + 'Output.txt', 'w') as output:
        output.write({'vf':vf,'rmin':rmin,'ert':ert,'vol':vh,'obj':oh})
    fig = plt.figure()
    fig.add_subplot(1, 2, 1).plot(vh, color='tab:blue')
    fig.add_subplot(1, 2, 2).plot(oh, color='tab:red')
    mddb.saveAs('Output_files\\' + 'Final_design')
