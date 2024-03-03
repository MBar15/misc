# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:51:13 2023

@author: Sotol
"""
import numpy as np
import pandas as pd


def ReadBulkDAT(name="bulk.dat"):
    f = open(name,"r")    
    content = (f.readlines())    
    GRID = []
    EleNoTri = []
    SPC = []
    Force = []
    for i in range(len(content)):
        CheckNOD = content[i].find("GRID")
        CheckTri = content[i].find("CPLSTS3")
        CheckSPC = content[i].find("SPC")
        CheckForce = content[i].find("FORCE")
        if CheckNOD != -1:
            GRID.append(content[i])
        if CheckTri != -1:
            EleNoTri.append(content[i])
        if CheckSPC != -1:
            SPC.append(content[i]) 
        if CheckForce != -1:
            Force.append(content[i])
    GRIDtable=pd.DataFrame([s.split() for s in GRID])    
    Nod = np.float64(GRIDtable.to_numpy()[:,[3,4,5]])
    if np.sum(Nod[:,2]) < 10**-2:
        Nod = Nod[:,[0,1]]
    Eletable=pd.DataFrame([s.split() for s in EleNoTri])    
    EleNo = np.intc(Eletable.to_numpy()[:,3:])
    
    SPCtable= pd.DataFrame([s.split() for s in SPC]) 
    SPCtable.dropna(inplace=True)
    FixNod = np.intc(SPCtable.to_numpy()[:,2])
    
    Forcetable= pd.DataFrame([s.split() for s in Force])
    Forcetable.dropna(inplace=True)
    
    ForceNod = np.intc(Forcetable.to_numpy()[:,2])

    return Nod,EleNo,FixNod,ForceNod