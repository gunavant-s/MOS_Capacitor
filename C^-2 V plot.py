# The acronym MOS stands for Metal oxide semiconductor. 
# An MOS capacitor is made of a semiconductor body or substrate, an insulator and a metal electrode called a gate.
# Practically the metal is a heavily doped n+ poly-silicon layer which behaves as a metal layer. 
# The dielectric material used between the capacitor plates is silicon dioxide (SiO2).

import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import sympy as sp
import math
from scipy import stats


# Function for finding the nearest number in a list
def closest(arr, K):

    res = arr[0]
    N = len(arr)
        # Traverse the array
    for i in range(1, N, 1):

            # If absolute difference
            # of K and res exceeds
            # absolute difference of K
            # and current element
        if (abs(K - res) > abs(K - arr[i])):
            res = arr[i]

        # Return the closest
        # array element
    return res

#file = "CV Data.csv"
file = input("Enter file name : ")

cnt =0

c = []
v = []

row = int(input("Enter row of data in file : "))
v_col = int(input("Enter col of V in file : "))
c_col = int(input("Enter col of C in file : "))

Cox = 1.7976931348623157e-308

with open(file,"r") as f:
    data = csv.reader(f)
    for i in data:
        cnt += 1
        if (cnt >= row):
            v.append(float(i[v_col - 1]))
            c.append(1/(float(i[c_col - 1])**2))
            Cox = max(Cox,float(i[c_col - 1]))

plt.plot(v,c,'-o')
plt.xlabel("Voltage")
plt.ylabel("1/C2")
plt.show()


start_idx = 0
end_idx = 0

while True:
    #print("-2.0 to -1.4 best for straight line")
    start_val = float(input("Enter the starting val of voltage : "))
    end_val = float(input("Enter the ending val of voltage : "))
    start_idx = v.index(closest(v,start_val)) - 3
    end_idx = v.index(closest(v,end_val)) - 3
    
    # -2.0 to -1.4 best for straight line
    
    # start_idx =  v.index(min(v, key = lambda x1 : abs(start_val - x1))) - 3
    # end_idx =  v.index(max(v, key = lambda x1 : abs(start_val - x1))) - 3
    print(f"Starting idx : {start_idx}\nEnding index : {end_idx}\n")
    plt.plot(v[start_idx:end_idx],c[start_idx:end_idx],'-o')
    plt.show()
    again = int(input("Press 1 to try again else 0 : "))
    if again == 1:
        pass
    else:
        break

curve = np.polyfit(v[start_idx:end_idx],c[start_idx:end_idx],1)
poly = np.poly1d(curve)
temp = [] 

for i in v[start_idx:end_idx]:
    temp.append(poly(i))
    
plt.plot(v[start_idx:end_idx],c[start_idx:end_idx],'-o',color='r')
plt.plot(v[start_idx:end_idx],temp[:],'--',color='b')
plt.legend(["Actual Data","Curve Fit"])

xs = sp.symbols('x')
f_lst = 0

for i in range(len(curve)):
    f_lst += curve[i]*pow(xs,(len(curve)-1-i))

df = sp.diff(f_lst, xs)
slope = df.subs(xs, 1).evalf() #SINCE LINEAR
print(slope)
# print((c[start_idx + 4] - c[start_idx])/(v[start_idx + 4] - v[start_idx]))

k = 1.380*math.pow(10,-23) # Boltzman const
q = 1.6*math.pow(10,-19) #charge of e
area = float(input("Enter area: "))
freq = float(input("Enter frequency(in Hz): "))
w = 2 * math.pi * freq
a_sq = area**2
temperature = float(input("Enter temperature : "))

Ԑo = 8.85418*math.pow(10,-14)
Ԑr = float(input("Enter Ԑr: "))
Ԑs = Ԑo*Ԑr

N = 2/(q*Ԑs*a_sq*slope)
LD = math.sqrt((Ԑs*k*temperature)/((q**2)*abs(N)))

Cs = Ԑs * area / LD
Cfb = (Cox*Cs)/(Cox + Cs)

Vfb = v[c.index(closest(c,Cfb))]

Ⴔb = -(k*temperature/q)*(math.log(abs(N)/9650000000))

Ⴔm = float(input("Enter Ⴔm : "))
ꭕ = float(input("Enter ꭕ : "))
Eg = float(input("Enter Eg : "))


#n-type for given excel sheet
print()

if(input("Enter which type (n or p) : ").lower() == "n"):
    Ⴔms = Ⴔm - (ꭕ + (Eg/2) - Ⴔb) # n-type
else:
    Ⴔms = Ⴔm - (ꭕ + (Eg/2) + Ⴔb) # p-type


Qeff = Cox*(Ⴔms - Vfb)/(q*area)
Neff = Qeff/q

print("\nSlope : ",slope)
print("N : ",N)
print("|N| : ",abs(N))
print("LD : ",LD)
print("Cox : ",Cox)
print("Cs : ",Cs)
print("Cfb : ",Cfb)
print("Vfb : ",Vfb)
print("Ⴔb : ",'%e' % Ⴔb)
print("Ⴔm : ",Ⴔm)
print("Ⴔms : ",'%e' % Ⴔms)
print("ꭕ : ",ꭕ)
print("Eg : ",Eg)
print("Qeff : ",'%e' % Qeff)
print("Neff : ",Neff)
print()


constants = dict()
constants[1] = area
constants[2] = freq
constants[3] = temperature
constants[4] = Ԑr
constants[5] = Cox
constants[6] = Ⴔm 
constants[7] = ꭕ
constants[8] = Eg

choice_dict = dict()

choice_dict[1] = "area"
choice_dict[2] = "freq"
choice_dict[3] = "temperature"
choice_dict[4] = "Ԑr"
choice_dict[5] = "Cox"
choice_dict[6] = "Ⴔm "
choice_dict[7] = "ꭕ"
choice_dict[8] = "Eg"

while True:
    again = int(input("Enter 1 to change constants else 0 : "))
    if again == 1:
        for k,d in choice_dict.items():
            print(f"{k} : {d}")
        while True:
            choice = int(input("Enter choice of constant to change : "))
            constants[choice] = float(input(f"Enter the value of {choice_dict[choice]} : "))
            # globals().get(choice_dict[choice], '') = constants[choice]
            
            if choice == 1:
                area = constants[choice]
            elif choice == 2:
                freq = constants[choice]
            elif choice == 3:
                temperature = constants[choice]
            elif choice == 4:
                Ԑr = constants[choice]
            elif choice == 5:
                Cox = constants[choice]
            elif choice == 6:
                Ⴔm = constants[choice]
            elif choice == 7:
                ꭕ = constants[choice]
            elif choice == 8:
                Eg = constants[choice]
            else:
                print("invalid")
            
            others = int(input("Any other constants to change 1 if yes else 0 : "))
            if others == 1:
                pass
            else:
                break
                
        w = 2 * math.pi * freq
        a_sq = area**2
        Ԑs = Ԑo*Ԑr

        N = 2/(q*Ԑs*a_sq*slope)
        LD = math.sqrt((Ԑs*k*temperature)/((q**2)*abs(N)))

        Cs = Ԑs * area / LD
        Cfb = (Cox*Cs)/(Cox + Cs)
        Vfb = v[c.index(closest(c,Cfb))]
        Ⴔb = -(k*temperature/q)*(math.log(abs(N)/9650000000))

        print()
        
        if(input("Enter which type (n or p) : ").lower() == "n"):
            Ⴔms = Ⴔm - (ꭕ + (Eg/2) - Ⴔb) # n-type
        else:
            Ⴔms = Ⴔm - (ꭕ + (Eg/2) + Ⴔb) # p-type


        Qeff = Cox*(Ⴔms - Vfb)/(q*area)
        Neff = Qeff/q
        print("\nSlope : ",slope)
        print("N : ",N)
        print("|N| : ",abs(N))
        print("LD : ",LD)
        print("Cox : ",Cox)
        print("Cs : ",Cs)
        print("Cfb : ",Cfb)
        print("Vfb : ",Vfb)
        print("Ⴔb : ",'%e' % Ⴔb)
        print("Ⴔm : ",Ⴔm)
        print("Ⴔms : ",'%e' % Ⴔms)
        print("ꭕ : ",ꭕ)
        print("Eg : ",Eg)
        print("Qeff : ",'%e' % Qeff)
        print("Neff : ",Neff)
        
        
    else:
        break

Gm = 1
Cm = 1
Gp_w_max = (w*Gm*Cox*Cox)/((Gm**2) + (w*w*(Cox - Cm)*(Cox - Cm)))
#Dit = (2 / LD)*((area*area/Cox)/((area*area/Cs)**2+(1-slope/area)**2))
Dit = 2.5*Gp_w_max/(q*area)
denom = ((Gm/(w*Cox))**2)+(1-Cm/Cox)**2
DIT = Dit/denom
