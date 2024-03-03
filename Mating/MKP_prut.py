import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

def Setup(ALFA=0):
    """Funkce načítající sít a formátující dané data pro řešič"""

    # Vstupní parametry úlohy, mat zatížení atd. 
    Young = 210.0e3
    Poisson = 0.3
    S = 1000.0
    SilaVelikost = 10.0e3
    
    # uzly 
    Nodals = np.array([
                        [0.0,0.0],
                        [1000.0,0.0],
                        [1000.0,1000.0]
                ])
        
    # elmenty (jaké uzly tvoři element) 
    EleNo = np.array([
                        [1,2],
                        [2,3],
                        [1,3]
                ])
    


    ET = EleNo.flatten() # převádí matici na pole ---------------------------------------
    # stupne volnosti elementu (levy xy a pravy xy kazdeho elmentu) ---------------------------------------
    DOF = np.concatenate([[2*ET-1, 2*ET]], axis=0).T #---------------------------------------
    DOF = DOF.reshape([np.size(EleNo,0), 2*np.size(EleNo,1)]) - 1 #---------------------------------------
    print(DOF)
    DOFu = np.unique(DOF.flatten()) # globalni stupne volnosti 6 (2*3 elementy) #---------------------------------------
    nDOF = np.max(DOFu)+1           # Counter stupnu volnosti  jedno vetsi #---------------------------------------
    FixDOF = np.array([0,1,3])      # stupne volnosti ktere jsou zamknute (locked) O.P. #---------------------------------------
    FreeDOF = np.setdiff1d(DOFu, FixDOF)    # STUPNE VOLNOSTI - O.P #---------------------------------------
    F = np.zeros([nDOF,1])                  # pole sil #---------------------------------------
    ForceDOF = [4,5] #np.array([])          # stupne volnosti sil  #---------------------------------------
    F[ForceDOF] = SilaVelikost  * np.array([np.cos(ALFA),np.sin(ALFA)]).reshape(2,1)    # pythagorova veta sila na stupne volnosti

        
    return {'Young':Young,'Poisson':Poisson,'Nodals':Nodals,'EleNo':EleNo,'DOF':DOF,'nDOF':nDOF,'FreeDOF':FreeDOF,'Sila':F,'Prurez':S}



def StiffnessBAR(Nodal,E,Area): # Uzly prvku, Younguv modul, plocha prurezu
    'Matice tuhosti elemetu lokalni a transformace do globalni SS'
    # velikost elementu
    A = Nodal[0,:]
    B = Nodal[1,:]
    V = B-A 
    Delka = np.linalg.norm(V) # pythagorova veta
    #Delka2 = (V.T @ V) **(1/2)


    c = V[0] / Delka    # COS
    s = V[1] / Delka    # SIN
    
    # tranformacni cleny do transformacni matice
    cc = c*c
    ss= s*s
    cs = c*s
    
    k = E*Area/Delka * np.array([[cc,cs,-cc,-cs],
                             [cs,ss,-cs,-ss],
                             [-cc,cs,cc,cs],
                             [-cs,-ss,cs,ss]
                             ])
    
    return k

def GlobalStiffness(Nodals,EleNo,DOF,nDOF,E,Prurez):
    """Globalni matice tuhosti"""

    K = np.zeros((nDOF,nDOF))   # Prazdna matice tuhosti #---------------------------------------
    LokalizacniTabulka = np.zeros((0,3)) #---------------------------------------
    for i in np.arange(0,np.size(EleNo,0)):  #---------------------------------------
        ke = StiffnessBAR(Nodals[EleNo[i,:]-1,:],E,Prurez) #---------------------------------------
# =============================================================================
#         rw = 0
#         for j in DOF[i,:]:
#             cl = 0
#             for k in DOF[i,:]:
#                 K[j,k] = K[j,k] + ke[cl,rw]
#                 cl = cl+1
#             rw=rw+1
# =============================================================================
        Ke = ke.flatten() #---------------------------------------
        cl = np.kron(np.ones([np.size(DOF[i,:])]), DOF[i,:]) #---------------------------------------
        rw = np.kron(DOF[i,:],np.ones([np.size(DOF[i,:])])) #---------------------------------------
        RCK = np.stack([rw,cl,Ke],axis=0).T #---------------------------------------
        LokalizacniTabulka = np.concatenate((LokalizacniTabulka,RCK), axis=0) #---------------------------------------
        
    K = sps.coo_matrix((LokalizacniTabulka[:,2],(LokalizacniTabulka[:,0],LokalizacniTabulka[:,1])),shape=(nDOF,nDOF)) #---------------------------------------
        
    return (K+K.T)/2 #---------------------------------------

def main(ALFA):

    # Read vstupni parametry
    SET = Setup(ALFA) #---------------------------------------
    E = SET['Young'] #---------------------------------------
    Prurez = SET['Prurez'] 
    Nodals = SET['Nodals'] #---------------------------------------
    EleNo = SET['EleNo'] #---------------------------------------
    DOF = SET['DOF'] #---------------------------------------
    nDOF = SET['nDOF'] #---------------------------------------
    freeDOF = SET['FreeDOF'] #---------------------------------------
    F = SET['Sila'] #---------------------------------------
    
    K = GlobalStiffness(Nodals,EleNo,DOF,nDOF,E,Prurez) # sestavenÁ GLOBÁLNÍ MATICE TUHOSTI
    K0 =K #---------------------------------------
    K = K.toarray() #---------------------------------------
    # APLIKACE OKRAJOVYCH PODMINEK
    Krow = K[freeDOF,:]  #---------------------------------------
    Kfree = Krow[:,freeDOF] #---------------------------------------

    u = np.zeros((nDOF,1))  # pole posuvu #---------------------------------------

    # reseni zakladni rovnice
    invK = np.linalg.inv(Kfree) #---------------------------------------
    u[freeDOF] = invK @ F[freeDOF] #---------------------------------------

    Ksparse = K0.tocsr()[freeDOF,:].tocsr()[:,freeDOF]
    usparse = sps.linalg.spsolve(Ksparse, F[freeDOF])

    uNodal = u.reshape((np.size(Nodals,0),2)) #---------------------------------------
    scaleValue = 500 
    plt.quiver(Nodals[:,0],Nodals[:,1], scaleValue*uNodal[:,0],scaleValue*uNodal[:,1],scale_units='xy',scale=1)
    plt.xlim([-100,1200])
    plt.ylim([-100,1200])
    plt.grid(0.95)
    
    #plt.figure()
    plt.scatter(Nodals[:,0], Nodals[:,1],alpha=0.5,s=50,c='k')
    barva = np.linalg.norm(uNodal,axis=1)
    plt.scatter(Nodals[:,0]+scaleValue*uNodal[:,0], Nodals[:,1]+scaleValue*uNodal[:,1],s=50,c=barva)
    plt.show()

    Sigma = napeti(Nodals,EleNo,DOF,E,u) #---------------------------------------
    #print(Sigma)
    
    
    return np.max(abs(Sigma)) #---------------------------------------


def napeti(Nodals,EleNo,DOF,E,u):
    Sigma = np.zeros((np.size(EleNo,0),1))
    
    for i in np.arange(0,np.size(EleNo,0)):
        uEGlobal = u[DOF[i,:]]
        B, cA, sA = Bbar(Nodals[EleNo[i,:]-1,:])
    
        TM = np.array([[cA,sA,0,0],[0,0,cA,sA]])
        uE = TM @ uEGlobal
        #print(E* (B @ uE))
        Sigma[i] = E* (B @ uE)
    return Sigma

def Bbar(Nodal):
    A = Nodal[0,:]
    B = Nodal[1,:]
    V = B-A
    Delka = np.linalg.norm(V)
    
    cA = V[0]/Delka
    sA = V[1]/Delka
    
    return np.array([-1/Delka,1/Delka]), cA, sA

if __name__ == '__main__':
    MAXSigma = 0
    #for uhel in np.arange(0,360,30):
    ALFA = np.radians(30)
    Sigma = main(ALFA)
        
        #if Sigma > MAXSigma:
        #    uhelMax = uhel
        #    MAXSigma = Sigma
    
    # print(f'Maximalni napeti {MAXSigma:.2f} je pri uhlu {uhelMax:.1f} stupnu')