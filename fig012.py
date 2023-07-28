import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
import time, traceback
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']

rg = 200
m0 = 0.1
rho = 1000
rhoa = 1.2250
S0 = 4.792e-4
S = 59.447e-4
H = 0.25233
ht0 = 0
p0 = ATM = 1.01325 * 100000
pt0 = 3 * ATM
g = 9.81
cd = 0.4
gamma = 1.4

'''
m0 = 0.1
rho = 1000
rhoa = 1.2250
S0 = 0.0006157
S = 0.007162
H = 0.33
ht0 = 0.33 / 3
p0 = ATM = 1.01325 * 100000
pt0 = 3 * ATM
g = 9.81
cd = 0.4
'''


def fcondi(argp):
    return((1-(10e-3))*(argp/ATM-1)/(-9)-(10e-3))

condi = fcondi(pt0)

def fp(argh):
    return(pt0 *(((H-ht0) / (H-argh))**gamma))

def fm(argh):
    return(m0 + rho * S * argh)


def vabsv(argv):
    return (argv / abs(argv)) if argv != 0 else 1


newv=[]
dydt = [0,0,0,0]
lastt = 0
lastv = 0
oof = 1e5
a = 0
def fun(t, y):
    global lastt,lastv,a,oof
    #print(t,y)
    h = y[0]
    v = y[1]
    #print(t)
    #input()
    #newv.append(v)
    pos_y = y[2]
    U = y[3]
    #print(U)
    p = fp(h)
    m = fm(h)
    A = (2 * S0) / ( m * (1-(S0**2/S**2)))
    '''
    dydt[1] = (A * (p + rho * g * h - p0) - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m) / (1-A*rho*h)
    a = dydt[1]
    dydt[0] = - (S0 / S) * np.sqrt(abs((2*(p+rho*(g+a)*h-p0))/(rho*(1-(S0**2/S**2)))))
    dydt[2] = v
    '''
    if t < oof:
        if h > 0:
            if pos_y <= 0:
                if (A * (p + rho * g * h - p0) - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m) / (1-A*rho*h) < 0:
                    dydt[1] = 0
                    a = 0
                else:
                    dydt[1] = (A * (p + rho * g * h - p0) - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m) / (1-A*rho*h)
                    a = dydt[1]
                if v < 0:
                    dydt[2] = 0
                else:
                    dydt[2] = v
            else:
                dydt[1] = (A * (p + rho * g * h - p0) - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m) / (1-A*rho*h)
                a = dydt[1]
                dydt[2] = v
            dydt[0] = - (S0 / S) * np.sqrt(abs((2*(p+rho*(g+a)*h-p0))/(rho*(1-(S0**2/S**2)))))
        else:
            dydt[1] = - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m
            a = dydt[1]
            dydt[2] = v
            dydt[0] = 0
        if dydt[0] >= condi:
            #pass
            oof = t
    else:
        if h > 0:
            if pos_y <= 0:
                if - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m < 0:
                    dydt[1] = 0
                    a = 0
                else:
                    dydt[1] = - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m
                    a = dydt[1]
                if v < 0:
                    dydt[2] = 0
                else:
                    dydt[2] = v
            else:
                dydt[1] = - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m
                a = dydt[1]
                dydt[2] = v
            dydt[0] = 0
        else:
            dydt[1] = - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m
            a = dydt[1]
            dydt[2] = v
            dydt[0] = 0
    #dydt[2] = v if not(pos_y <= 0 and dydt[1] < 0) else 0
    
    '''
    #u = np.sqrt((2*(p + rho * (g+a) * h - p0) / rho)-v**2)
    # dh/dt = [-((S0/S) * np.sqrt(abs((fp(h)+ rho * g * h - p0) / rho)))]
    # dv/dt = (p + rho * g * h - p0) * S0 / m - g - 0.5 * cd * rhoa * S * v * v
    # dydt = [-((S0/S) * np.sqrt(abs((p + rho * g * h - p0) / rho))), ((p + rho * g * h - p0) * S0 / m - g - 0.5 * cd * rhoa * S * v * v / m), v]
    if h <= 0:
        #dydt[0] = 0
        #dydt[1] = - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m
    else:
        #dydt[0] = -((S0/S)*(p + rho * (g+a) * h - p0)/(rho*u))
        #dydt[1] = ((rho * (u**2) * S0 / m) - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m)  
    u = np.sqrt((2*(p + rho * (g+a) * h - p0) / ((1-(S0**2/s**2))*rho))
    if pos_y <= 0:
        if dydt[1] < 0:
            dydt[1] = 0
        if dydt[2] < 0:
            dydt[2] = 0
    dydt[2] = v if not(pos_y <= 0 and dydt[1] < 0) else 0
    '''
    #print(a)
    '''
    print("-----------------------------")
    print("h = ", h)
    print("v = ", v)
    print("t = ", t)
    print("m = ", m)
    print("p = ", p)
    print('u = ', u)
    print('y = ', pos_y)
    print("dh/dt = ",dydt[0])
    print("dv/dt = ",dydt[1])
    print("dy/dt = ",dydt[2])
    
    print("S0/S=", S0/S)
    print("A=", np.sqrt(abs((fp(h)+ rho * g * h - p0) / rho)))
    print("B=", (p+ rho * g * h - p0))
    print("p()=", p)
    print("rho*g*h=", rho * g * h)
    '''
    #input()
    
    return dydt


fig0 = plt.figure('fig0',figsize=(4,2), dpi=150)
ax0 = fig0.add_subplot(111)
plt.xlabel("时间,s")
plt.ylabel("水高,m")
plt.axis([0,0.25,0,H])
plt.grid()
plt.tight_layout()

fig1 = plt.figure('fig1',figsize=(4,2), dpi=150)
ax1 = fig1.add_subplot(111)
plt.xlabel("时间,s")
plt.ylabel("速度,m/s")
plt.axis([0,2.5,-5,30])
plt.grid()
plt.tight_layout()

fig2 = plt.figure('fig2',figsize=(4,2), dpi=150)
ax2 = fig2.add_subplot(111)
plt.xlabel("时间,s")
plt.ylabel("高度,m")
plt.axis([0,5,0,30])
plt.grid()
plt.tight_layout()


for i in [0.3,0.4,0.5,0.6]:
    try:
        dydt = [0,0,0,0]
        a=0
        oof = 1e5
        lastt = 0
        lastv = 0
        ht0 = H * i
        # 初始条件
        y0 = [ht0, 0, 0,0]
        #print(y0)
        rge = 10
        yy = solve_ivp(fun, (0, rge), y0, method='Radau',t_eval = np.arange(0,rge,0.01) )
        #xx = solve_ivp(fun, (0, rg), y0, method='Radau')

        t = yy.t
        data = yy.y
        max_y = np.amax(data[2])
        where_y = np.argmax(data[2])
        print(t[where_y], max_y)
        #plt.plot(t, data[0, :],color='blue')
        #plt.plot(t, data[1, :],color='red')
        ax0.plot(t, data[0, :],color='black',linewidth=1)
        ax1.plot(t, data[1, :],color='black',linewidth=1)
        ax2.plot(t, data[2, :],color='black',linewidth=1)
        #t2 = xx.t
        #data2 = xx.y
    except:
        traceback.print_exc()

plt.show()
