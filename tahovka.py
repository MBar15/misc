import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# sample info 
l1 = 10
d1 = 5

l2 = 10
d2 = 5

l3 = 10
d3 = 5

S = np.pi*d1**2/4 

# Reading file 
file = np.loadtxt('C:/Users/Michal/Desktop/Magisterský/MatNaPoIngUloh/TensileTest1_MS300_10mmPerMin_D200422_T103312.csv',float,delimiter=';',usecols=(1,3))
# Reading proper data
Force, elongation = file[:,0],file[:,1]

# Deleting data wrong data
    # Negative force
index = []
for i in range(len(Force)):
    if Force[i] < 0:
        index.append(i)

Force = np.delete(Force,index)
elongation = np.delete(elongation,index)

    # data after tearing
stop_elongation = 0
e_prev = 0
for i in range(elongation.size):
    if elongation[i]-e_prev < 0 and i> 3/4 * elongation.size:
        #print(e-e_prev , e_prev, e)
        stop_elongation = i
        break
    e_prev = elongation[i]

delete_index = [i for i in range(stop_elongation,elongation.size)]

Force = np.delete(Force,delete_index)
elongation = np.delete(elongation,delete_index)

# Calculating stress and strain
Stress = np.ndarray(Force.size)
Strain = np.ndarray(Force.size)
ind = 0
for i in range(Force.size):
    Stress[ind] = Force[i]/S
    Strain[ind] = elongation[i]/l1
    ind += 1



# young module
# find linear part
E_final = 0
q_final = 0
R_final = 0
max_data = 10
step = 100
while True:
    del_dat = [x for x in range(max_data,Stress.size)]
    newStress = np.delete(Stress,del_dat)
    newStrain = np.delete(Strain,del_dat)
    x = np.array(newStrain).reshape(-1,1)
    y = np.array(newStress)
    regression_model = LinearRegression()
    regression_model.fit(x,y)


    E_final = regression_model.coef_
    q_final = regression_model.intercept_
    R_final = regression_model.score(x,y)
    #print(regression_model.score(x,y), max_data)
    if regression_model.score(x,y) > 0.9998 : break

    max_data += step


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Tahová zkouška výsledky')

ax1.plot(Strain,Stress,'-')
ax1.set_ylabel('Stress [MPa]')
ax1.set_xlabel('Strain [-]')

ax2.plot(x,y,
          x,E_final*x+q_final)
ax2.set_ylabel('Stress [MPa]')
ax2.set_xlabel('Strain [-]')
ax2.text(0.00015, 400, f'R^2: {np.round(R_final,5)}', fontsize=10)
ax2.text(0.00015, 380, f'E: {np.round(E_final[0],2)}', fontsize=10)
ax2.text(0.00015, 360, f'y={np.round(E_final[0],2)}*x + {np.round(q_final,2)}', fontsize=10)
print(type(E_final))
plt.show()

print(R_final)