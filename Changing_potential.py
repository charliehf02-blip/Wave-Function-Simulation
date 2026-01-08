
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import simps


N=1000


gamma_sqrd = 1000

l = 1/N
psi=np.zeros(N)
psi[0]=0 
psi[1] = 1e-4





potential_arr = []

x_bar = np.arange(0, 1 , 0.001)
for i in range(len(x_bar)):
    potential_arr.append(8*(x_bar[i]-0.5)**2 - 1)

#plt.plot(x_bar,potential_arr)



def psi_array(eps):
    k_sqrd =np.zeros(N)
    for i in range(len(potential_arr)):
        k_sqrd[i] = gamma_sqrd * ( eps - potential_arr[i]) 
    n=1
    while (n < N-1):
        psi[n+1] = (2*(1-5/12 * l**2 * k_sqrd[n])*psi[n] - (1 + 1/12 * l**2 * k_sqrd[n-1])*psi[n-1])/(1+1/12*l**2*k_sqrd[n+1])
        #psi[n+1] = psi[n] + 1
        n+=1
    return psi




def psi_norm(eps):
    psi_norm = psi_array(eps)/(np.sqrt(simps(psi**2,dx=l)))
    return psi_norm




#plt.plot(psi_array(-0.9105468750000001))
#plt.plot(np.zeros(N),'r--' , label='zeros')





def find_eigen(eps,tol):

    delta_eps = 0.001 
    while abs(delta_eps)>tol:
        sign_change = psi_array(eps)[N-1] * psi_array(eps+ delta_eps)[N-1]
        eps+=delta_eps
        if sign_change<0:  
            #print("less than 0: ", sign_change, "delta_eps is: ", delta_eps, "epsilon is: ", eps)
            delta_eps = -delta_eps / 2
    return eps

tol = 0.000001

lst = np.arange(-1,4,0.1)
temp=[]
for i in range(len(lst)):
    #print(find_eigen(lst[i], tol))
    temp.append(find_eigen(lst[i], tol))


eigenvalues =[]
for i in range(len(temp)-1):
    if abs(temp[i+1] - temp[i]) >0.1:
        eigenvalues.append(temp[i])
        
print("eigenvalyes are" , eigenvalues) 
print("amoubnt of " , len(eigenvalues))       
#for i in range(len(eigenvalues)-1):
#    plt.plot(psi_norm(eigenvalues[i+1]))
    
    
del_eig=[]

for i in range(len(eigenvalues)-1):
    del_eig.append(abs(eigenvalues[i+1] - eigenvalues)[i])






plt.plot(del_eig)
plt.title('difference between adjacent eigenvalues')
plt.xlabel("n")
plt.ylabel("delta e")
plt.yscale('log')
plt.show()



firs_ten_eig = np.zeros(10)
for i in range(10):
    firs_ten_eig[i] = eigenvalues[i]

plt.plot(firs_ten_eig)
plt.title('First Ten Eigenvalues')
plt.xlabel("n")
plt.ylabel("e")

plt.show()


eq = np.zeros(N)
plt.plot(psi_norm(eigenvalues[12]))
plt.plot(eq, 'r--')





#print(v_pot)

"""

eps_tries = np.arange(-1,4,0.1)
#print(v_pot)

temp =[]

for i in range(len(eps_tries)):
    temp.append(find_eigen(,tol))
    
eigenvalues =[]
for i in range(len(temp)-1):
    if abs(temp[i+1] - temp[i]) >0.1:
        eigenvalues.append(temp[i])
        

print("eigenvalues" , eigenvalues)

 


plt.title("Wave Function in Harmonic potential")
plt.xlabel("N-iterations")
plt.ylabel("Psi value")
plt.show()


"""



def psi_pp(eps):
    psi = psi_norm(eps)
    
    psi_pp = np.zeros(N)
    for n in range(N-1):
        if n !=0:
            if n != N-1:
                psi_pp[n] = (psi[n-1] - 2*psi[n] + psi[n+1]) /(l**2)
        
    return psi_pp



plt.show()
plt.clf()

#plt.plot(psi_pp_analytic(eigenvalues[0]))


def delp(eps):
    delp = np.sqrt( -1*simps(psi_norm(eps)*psi_pp(eps), dx = l))
    return delp
   
x = np.arange(0,1,l)

def delx(eps):
    delx = np.sqrt( simps(x**2 * psi_norm(eps)**2,dx=l) - (simps(x*psi_norm(eps)**2,dx=l)**2))
    return delx

print("delp is " , delp(eigenvalues[0]))
print("delx is " , delx(eigenvalues[0]))
print("delxdelp is " , delx(eigenvalues[0])*delp(eigenvalues[0]))



delx_arr=[]
delp_arr=[]
delxp_arr =[]

for i in range(len(eigenvalues)):
    print("delxdelp is " , delx(eigenvalues[i])*delp(eigenvalues[i]))
    delx_arr.append(delx(eigenvalues[i]))
    delp_arr.append(delp(eigenvalues[i]))
    delxp_arr.append( delx(eigenvalues[i])*delp(eigenvalues[i]))


plt.plot(delx_arr, 'r--', label='delta x')
plt.plot(delp_arr,  'g--', label='delta p')
plt.plot(delxp_arr,  'b', label='del x*del p')
    


plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


plt.show()












