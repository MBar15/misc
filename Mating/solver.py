import numpy as  np
import scipy.sparse as sps

def StiffnessCST(Nodal,E,nu,t):
    D = E/(1-nu**2) * np.array([
                                [1,nu,0],
                                [nu,1,0],
                                [0,0,(1-nu)/2]
    ])

    Node1 = Nodal[0,:]
    Node2 = Nodal[1,:]
    Node3 = Nodal[2,:]
    b1 = Node2[1]-Node1[1]
    b2 = Node3[1]-Node1[1]
    b3 = Node1[1]-Node2[1]
    c1 = Node3[0]-Node2[0]
    c2 = Node1[0]-Node3[0]
    c3 = Node2[0]-Node1[0]

    A = 1/2*np.linalg.det([
                        [1,Node1[0],Node1[1]],
                        [1,Node2[0],Node2[1]],
                        [1,Node3[0],Node3[1]]
    ])

    B = 1/(2*A) * np.array([
                            [b1,0,b2,0,b3,0],
                            [0,c1,0,c2,0,c3],
                            [c1,b1,c2,b2,c3,b3]
    ])

    k = A*t* np.transpose(B) @ D @ B
    return k

def GlobalStiffness(Nodals,EleNo,DOF,nDOF,E,nu,t):
    """Globalni matice tuhosti"""

    K = np.zeros((nDOF,nDOF))   # Prazdna matice tuhosti
    LokalizacniTabulka = np.zeros((0,3))
    for i in np.arange(0,np.size(EleNo,0)):
        ke = StiffnessCST(Nodals[EleNo[i,:]-1,:],E,nu,t)
        Ke = ke.flatten()
        cl = np.kron(np.ones([np.size(DOF[i,:])]), DOF[i,:])
        rw = np.kron(DOF[i,:],np.ones([np.size(DOF[i,:])]))
        RCK = np.stack([rw,cl,Ke],axis=0).T
        LokalizacniTabulka = np.concatenate((LokalizacniTabulka,RCK), axis=0)
        
    K = sps.coo_matrix((LokalizacniTabulka[:,2],(LokalizacniTabulka[:,0],LokalizacniTabulka[:,1])),shape=(nDOF,nDOF))
        
    return (K+K.T)/2


def Stress(Nodals,EleNo,DOF,E,u): #---------------------------------------
    Sigma = np.zeros((np.size(EleNo,0),1))
    
    for i in np.arange(0,np.size(EleNo,0)):
        uEGlobal = u[DOF[i,:]]
        #print(uEGlobal)
        D, B = CST(Nodals[EleNo[i,:]-1,:])
        uu = D @ B
        #print(uu)
        #print(uu @ uEGlobal)
        Sigma[i] =  uu @ uEGlobal
    return Sigma


def CST(Nodal):
    E = 2.1e6
    nu = 0.3
    D = E/(1-nu**2) * np.array([
                                [1,nu,0],
                                [nu,1,0],
                                [0,0,(1-nu)/2]
    ])

    Node1 = Nodal[0,:]
    Node2 = Nodal[1,:]
    Node3 = Nodal[2,:]
    b1 = Node2[1]-Node1[1]
    b2 = Node3[1]-Node1[1]
    b3 = Node1[1]-Node2[1]
    c1 = Node3[0]-Node2[0]
    c2 = Node1[0]-Node3[0]
    c3 = Node2[0]-Node1[0]

    A = 1/2*np.linalg.det([
                        [1,Node1[0],Node1[1]],
                        [1,Node2[0],Node2[1]],
                        [1,Node3[0],Node3[1]]
    ])

    B = 1/(2*A) * np.array([
                            [b1,0,b2,0,b3,0],
                            [0,c1,0,c2,0,c3],
                            [c1,b1,c2,b2,c3,b3]
    ])


    return D, B