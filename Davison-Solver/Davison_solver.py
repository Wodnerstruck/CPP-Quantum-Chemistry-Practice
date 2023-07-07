#Davison Diag
import numpy as np
import math
tol = 1e-8
import time
n = 1500
mmax = n // 2
sparsity = 0.001
A = np.zeros((n,n))
for i in range(0,n):
    A[i,i] = i + 1
A = A + sparsity*np.random.randn(1000,1000)
A = (A.T + A) / 2

k = 8 # number of initial guess vectors
eig = 4 # number of eigenvalues to solve
t = np.eye(n,k) # set k unit initial guess vectors
V = np.zeros((n,n))
I = np.eye(n)

#Begin Davidson
start_davidson = time.time()

for m in (k,mmax,k):
    if m <= k:
        for i in range(0,m):
            V[:,i] = t[:,i]/np.linalg.norm(t[:i]) # set initial guess vectors
        theta_old = 1 # arbitrary value
    elif m > k:
        theta_old = theta[:eig]
    V[:,:m], R = np.linalg.qr(V[:,:m]) # orthogonalize V
    H = np.dot(V[:,:m].T,np.dot(A,V[:,:m])) # form H matrix
    THETA, U = np.linalg.eig(H) # solve H matrix
    index = np.argsort(THETA) # sort eigenvalues index
    theta = THETA[index] # sort eigenvalues
    u = U[:,index] # sort eigenvectors
    
    for i in range(0, k):
        q = np.dot((H - theta[i]*I),np.dot(V[:,:m], u[:i])) / (theta[i] - H[i,i]) # form correction vector
        V[:,m + i] = q
    norm = np.linalg.norm(theta[:eig] - theta_old) # check convergence
    if norm < tol:
        break

end_davidson = time.time()

print("davidson = ", theta[:eig],";",
    end_davidson - start_davidson, "seconds")

start_numpy = time.time()

E,Vec = np.linalg.eig(A)
E = np.sort(E)

end_numpy = time.time()


print("numpy = ", E[:eig],";",
     end_numpy - start_numpy, "seconds") 
        
             
