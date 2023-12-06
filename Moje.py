import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt
import ReadData
from solver import *



def Setup():
    """Funkce načítající sít a formátující dané data pro řešič"""

    # Material constants
    t = 10 # thickness mm
    E = 2.1e5 # Young modul MPa
    nu = 0.3 # Poisson ratio
    Force = 10.e3#10.0e3 # Force magnitude applied
    
    Nodals, EleNo, FixDOF, ForceDOF = ReadData.ReadBulkDAT(name='bulk.dat')
    ET = EleNo.flatten() # převádí matici na pole

    # stupne volnosti elementu (levy xy a pravy xy kazdeho elmentu) 
    DOF = np.concatenate([[2*ET-1, 2*ET]], axis=0).T
    DOF = DOF.reshape([np.size(EleNo,0), 2*np.size(EleNo,1)]) - 1
 
    DOFu = np.unique(DOF.flatten()) # globalni stupne volnosti 6 (2*3 elementy)
    nDOF = np.max(DOFu)+1           # Counter stupnu volnosti  jedno vetsi
    # Fixxing DOF #---------------------------------------
    clamped_boundary = []
    for i in range(np.size(Nodals,0)):
        if Nodals[i][0] == 0:
            clamped_boundary.append(i)
    clamped_boundary = np.array(clamped_boundary)
    FixDOF = np.concatenate((clamped_boundary,clamped_boundary+1),axis=0)    # stupne volnosti ktere jsou zamknute (locked) O.P.
    FreeDOF = np.setdiff1d(DOFu, FixDOF)    # STUPNE VOLNOSTI - O.P


    F = np.zeros([nDOF,1])                  # pole sil
    Force_boundary = []
    for i in range(np.size(Nodals,0)):
        if Nodals[i][0] == 500:
            Force_boundary.append(i)
    ForceDOF = np.array(Force_boundary)+1   # stupne volnosti sil 
    Force_Per_Node = Force/ForceDOF.size
    F[ForceDOF] = -Force_Per_Node  * np.ones((ForceDOF.size,1))  
    Nodalss = Nodals[Force_boundary]
    plt.scatter(Nodalss[:,0], Nodalss[:,1],alpha=0.5,s=50,c='k')
    plt.show()
    
    return {'Young':E,'Poisson': nu,'Nodals':Nodals,'EleNo':EleNo,'DOF':DOF,'nDOF':nDOF,'FreeDOF':FreeDOF,'Sila':F,'Thickness':t}


def main():

    # Read vstupni parametry
    SET = Setup()
    E = SET['Young']
    nu = SET['Poisson'] #---------------------------------------
    t = SET['Thickness'] #---------------------------------------
    Nodals = SET['Nodals']
    EleNo = SET['EleNo']
    DOF = SET['DOF']
    nDOF = SET['nDOF']
    freeDOF = SET['FreeDOF']
    F = SET['Sila']
    
    K = GlobalStiffness(Nodals,EleNo,DOF,nDOF,E,nu,t) # sestavenÁ GLOBÁLNÍ MATICE TUHOSTI #---------------------------------------
    K0 =K
    K = K.toarray()

    # APLIKACE OKRAJOVYCH PODMINEK
    Krow = K[freeDOF,:] 
    Kfree = Krow[:,freeDOF]

    u = np.zeros((nDOF,1))  # pole posuvu

    # reseni zakladni rovnice
    invK = np.linalg.inv(Kfree)
    u[freeDOF] = invK @ F[freeDOF]

    uNodal = u.reshape((np.size(Nodals,0),2))

    uMagnitude = np.zeros((np.size(uNodal,0))) #---------------------------------------
    for i in range(np.size(uNodal,0)): #---------------------------------------
        uMagnitude[i] = np.sqrt(uNodal[i][0]**2 + uNodal[i][1]**2) #---------------------------------------

    scaleValue = 500
    plt.quiver(Nodals[:,0],Nodals[:,1], scaleValue*uNodal[:,0],scaleValue*uNodal[:,1],scale_units='xy',scale=0.5)
    plt.xlim([-20,600])
    plt.ylim([-20,125])
    plt.grid(0.95)
    
    #plt.figure()
    plt.scatter(Nodals[:,0], Nodals[:,1],alpha=0.5,s=50,c='k')
    barva = np.linalg.norm(uNodal,axis=1)
    plt.scatter(Nodals[:,0]+scaleValue*uNodal[:,0], Nodals[:,1]+scaleValue*uNodal[:,1],s=50,c=barva)
    plt.show()

    Sigma = Stress(Nodals,EleNo,DOF,E,u)
    
    
    return np.max(abs(Sigma))



if __name__ == '__main__':
    Sigma = main()
    print(Sigma)
        
