#!/usr/bin/python
# programmer : aber
# usage:
import sys
import time
import numpy as Np
import scipy as Sp
from numpy import linalg as Nla
from scipy import linalg as Sla
from scipy import stats 
import math

#def K_means(data):
    

def Matrix_dot(A,B): #Not Np.dot
    #print A.shape,B.shape
    if A.shape[1] == B.shape[0]:
        pass
    else:
        print "Matrix line-size are not equal !!!"
        print A.shape,A,B.shape,B
        sys.exit()
    x = Np.zeros([A.shape[0],B.shape[1]])
    AT = Sp.transpose(A)
    #print AT
    for i in range(AT.shape[1]):
        for j in range(B.shape[1]):
            x[i,j] = Np.dot(AT[:,i],B[:,j])
            #print x[i,j]
    return x

def PDF(point,MU,T): # generate pdf of a point in a serious of Gaussian with covanriance matrix T
    if point.shape[0] == MU.shape[0]:
        pass
    else:
        print "the dimonsions of data and Gaussian is different !!"
        sys.exit()
    rpdf = Np.zeros((MU.shape[1],point.shape[1]))
    #print MU.shape,'rpdf',rpdf
    _minused = Np.zeros((MU.shape[0],point.shape[1]))
    for m in range(MU.shape[1]):
        _invT = Nla.inv(Matrix_dot(Si[m,0],Sp.transpose(Si[m,0])))
        DET = Sla.det(_invT)
        _outerS = 1/((math.sqrt(DET))*((2*math.pi)**(MU.shape[0]/2.0)))
        _minused = point - MU[:,m].reshape(MU.shape[0],1)
        #print _minused
        _innerMa = Matrix_dot(Sp.transpose(_minused),_invT)
        _innerMb = Matrix_dot(_innerMa,_minused)
        _trace = 0.0
        for j in range(_innerMb.shape[0]):
            _trace += _innerMb[j,j]
        rpdf[m,0] = (-0.5)* _trace
    #print m,'rpdf__:   ',rpdf
    Gau_pdf = _outerS*(Np.exp(rpdf))
    for i in range(Gau_pdf.shape[0]):
        if Gau_pdf[i,0]  < 1.0e-80:
            Gau_pdf[i,0] = 1.0e-80
    #print 'PDF is ok',Gau_pdf
    return Gau_pdf

def Get_alpha_sigma_table(Gau_pdf,alpha_table):#alpha_table has k row, 1 column
    k = alpha_table.shape[0]
    n = Gau_pdf.shape[1]
    AS_table = Np.zeros((k+1,n))
    for i in range(k):
        AS_table[i,:] = alpha_table[i,:]*Gau_pdf[i,:]
        AS_table[k,:] += AS_table[i,:]
    #print 'AS_table is OK',AS_table
    return AS_table

def Yang_step(AS_table):#Have problem here !!!!!!
    k = AS_table.shape[0]
    n = AS_table.shape[1]
    pEM = Np.zeros((k,n))
    for i in range(k-1):
        #print AS_table[i,:],AS_table[k,:]
        pEM[i,:] = (AS_table[i,:]/AS_table[k-1,:])
        pEM[k-1,:] += pEM[i,:]*Np.log(pEM[i,:])
    Pit =  Np.zeros((k,n))
    for i in range(k-1):
        Pit[i,:] = pEM[i,:]+pEM[i,:]*(Np.log(pEM[i,:])-pEM[k-1,:])
        Pit[k-1,:] += Pit[i,:]
    #print 'Yang step is OK',Pit
    return Pit

def Ying_step(X,Pit,ps):#ps: penalty score
    k = Pit.shape[0]-1
    #update alpha
    _temp_alpha = Np.zeros((k,2))
    _temp_alpha[:,0] = (1-ps)*alpha[:,0]+ps*Pit[0:k,0]/Pit[k,0]
    alpha[:,0] = _temp_alpha[:,0]
    alpha[:,0].reshape(k,1)
    #update covariance
    n = mu.shape[0]
    M1 = Np.zeros((n,n))
    for j in range(n):
        M1[j,j] = 1
    for i in range(k):
        _mined = X[:,0]-mu[:,i].reshape(n,1)
        _mined_2 = Matrix_dot(Nla.inv(Matrix_dot(Si[i,0],Sp.transpose(Si[i,0]))),_mined)
        _mined_3 = Matrix_dot(_mined_2,Sp.transpose(_mined))
        #print 'Pit[i,:]',Pit#'mine-3',_mined_3,(_mined_3-M1),Pit[i,0]
        _Gi = Pit[i,0]*(_mined_3-M1)
        #print Gi
        _si = Matrix_dot(Si[i,0],(M1+ps*_Gi))
        #print'Si:', _si
        if Sla.det(_si)<0:
            _si = M1
        Si[i,0] = _si
    #update mu
    _temp_mu = Np.zeros((mu.shape))
    for i in range(k):
        #print mu[:,i],X.reshape(1,n),Pit[i,0]
        _temp_mu[:,i] = mu[:,i]+ps*Pit[i,0]*(X.reshape(1,n)-mu[:,i])
        mu[:,i] = Sp.transpose(_temp_mu[:,i])
    return

def YoN_discard(Th,ALPHA,MU,SI):
    k = ALPHA.shape[0]
    for i in range(k):
        if ALPHA[i,0] < Th:
            print 'Deleting : ',ALPHA[i,:]
            _alpha = Np.zeros((k,0))
            _alpha = ALPHA+ALPHA[i,0]/(k-1)
            _alpha = Np.delete(ALPHA,i,axis=0)
            _mu = Np.delete(MU,i,axis=1)
            _Si = Np.delete(SI,i,axis=0)
            break
        else:
            _alpha = ALPHA
            _mu = MU
            _Si = SI
    return _alpha, _mu, _Si

def Line2nparray(Line):
    Li = Line[:-1].split('\t')
    Lnp = Np.zeros((1,len(Li)))
    for i in range(len(Li)):
        Lnp[0,i] = float(Li[i])
    Np.random.shuffle(Lnp)
    return Lnp

def Get_variance(si):
    _invT = Np.copy(si)
    for m in range(si.shape[0]):
        _invT[m,0] = Matrix_dot(Si[m,0],Sp.transpose(Si[m,0]))
    return _invT

def Himoney(Data,ALPHA,MU,SI):
    #print Data.shape,MU.shape,ALPHA.shape,SI.shape
    sigma_H = 0.0
    _Hi = Np.zeros(Data.shape)
    for i in range(ALPHA.shape[0]):
        #print i,Data.shape, Np.array([[MU[0,i]]]).shape
        _all_pdf = PDF(Data,Np.array([[MU[0,i]]]),SI[i,0])
        #print 'pdf',_all_pdf
        _AS_table = Get_alpha_sigma_table(_all_pdf,Np.array([[ALPHA[i,0]]]))
        #print _AS_table
        _tem_hi = _AS_table[0,:]/_AS_table[1,:]
        #print _tem_hi
        _Hi = _tem_hi*Np.log(_AS_table[0,:])
    #print 'aaa',_Hi
    for each in _Hi:
        #print each
        sigma_H = sigma_H + each
        #print sigma_H
    return sigma_H

start = time.time()
#File = open(sys.argv[1])
te = Np.loadtxt(sys.argv[1])
B = Np.log2(te)
print max(B),min(B),B.shape
dis = max(B)-min(B)
mu = Np.array([[min(B),0.25*dis+min(B),0.5*dis+min(B),0.75*dis+min(B),max(B)]])
bT1 = Np.array([[1.0],[1.0]])#,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
COV = Np.array([[bT1],[bT1],[bT1],[bT1],[bT1]])

Si = Np.copy(COV)
for i in range(COV.shape[0]):
    Si[i,0] = COV[i,0]
alpha = Np.array([[0.2],[0.2],[0.2],[0.2],[0.2]])
Hi = Np.zeros((1,5))
sigma_hi = 0.0
PS = 0.015
n = 0
while n <50:
    Np.random.shuffle(B)
    #Last_alpha = Np.copy(alpha)
    #Last_sigma_hi = sigma_hi
    for i in range(B.shape[0]):
        #print B[i]
        rpdf = PDF(Np.array([[B[i]]]),mu,Si)
        AS = Get_alpha_sigma_table(rpdf,alpha)
        Pit = Yang_step(AS)
        Ying_step(Np.array([[B[i]]]),Pit,PS)
        [alpha,mu,Si]=YoN_discard(0.03,alpha,mu,Si)
        variance = Get_variance(Si)
        Last_alpha = Np.copy(alpha)
        Last_sigma_hi = sigma_hi
        if n > 1 and i in[200,400,600,800]:
            if Last_alpha.shape == alpha.shape:
                sigma_hi = Himoney(Np.array([B]),alpha,mu,Si)
                print 'delta: ',sigma_hi,Last_sigma_hi,(sigma_hi-Last_sigma_hi),0.002*abs(Last_sigma_hi)
                if (abs(sigma_hi-Last_sigma_hi) - 0.002*abs(Last_sigma_hi))< 0 :
                    n = 100
                    break
                PS = PS*0.1
                print 'New PS: ',PS
    n += 1
    print 'The ',n,' rounds:'
print 'Final alpha',alpha
print 'Final mu',mu
print 'Final var',variance
end = time.time()
print end-start
