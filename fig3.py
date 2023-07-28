import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.integrate import solve_ivp
import numpy as np
import time, traceback
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']

condi = -1e-1

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

def timeoutput(SecToConvert):
    RemainingSec = SecToConvert % (3600)
    HoursGet = SecToConvert // (3600)
    MinutesGet = RemainingSec // 60
    RemainingSec %= 60
    return(str(HoursGet)+'h '+str(MinutesGet)+'min '+str("%.d"%RemainingSec)+'s')

def timepredict(argtime):
    global N
    return N*argtime

def fcondi(argp):
    return((1-(10e-3))*(argp-1)/(-9)-(10e-3))

def fp(argh):
    return(pt0 *(((H-ht0) / (H-argh))**gamma))

def fm(argh):
    return(m0 + rho * S * argh)

def vabsv(argv):
    return (argv / abs(argv)) if argv != 0 else 1

dydt = [0,0,0]
oof = oof0=1e5
a = 0
def fun(t, y):
    global a,oof
    #print(t,y)
    h = y[0]
    v = y[1]
    #print(t)
    #input()
    #newv.append(v)
    pos_y = y[2]
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
                tempa = (A * (p + rho * g * h - p0) - g - 0.5 * cd * rhoa * S * v * v * vabsv(v)/ m) / (1-A*rho*h)
                if  tempa < 0:
                    dydt[1] = 0
                    a = 0
                else:
                    dydt[1] = tempa
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
        if dydt[0] >= condi and oof == oof0:
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

# 初始条件
y0 = [ht0, 0, 0]
#print(y0)
rge = 10

def find_root(data):
    global rge
    lastone = 0
    startime = 0
    for t in range(1,int(rge/0.01)):
        if lastone == 0 and data[t]>0:
            starttime = t
            break
        lastone = data[t]
    lastone = 0
    for t in range(1,int(rge/0.01)):
        if lastone > 0 and data[t] <= 0:
            return(t-startime)
            break
        lastone = data[t]

h_eval = np.arange(0,100,1)
p_eval = np.arange(3,3.1,1)
N = len(h_eval) * len(p_eval)
#h_best = 0 * p_eval
h_best = []
alt_best = []
t_best = []
h_now = len(h_eval)*[np.NaN]
posi = 0
lasttime = 0
for j in p_eval:
    condi = fcondi(j)
    best_time_dur = 0
    time_dur = 0
    best_h = 0
    best_alt = 0
    best_t = 0
    dontprint = 0
    yesdraw = 0
    pt0 = j * ATM
    cnt = -1
    for i in h_eval:
        cnt += 1
        N -= 1
        try:
            ht0 = H*i*0.01
            a=0
            oof = 1e5
            y0 = [ht0, 0, 0]
            dydt = [0,0,0]
            yy = solve_ivp(fun, (0, rge), y0, method='Radau',t_eval = np.arange(0,rge,0.01) )
            #xx = solve_ivp(fun, (0, rg), y0, method='Radau')

            t = yy.t
            data = yy.y
            max_y = np.amax(data[2])
            where_y = np.argmax(data[2])
            time_dur = find_root(data[2]) * 0.01
            if not dontprint:
                #print("%.4f" % ht0,'&',"%.2f" % float(ht0/H),'&',"%.4f" % max_y, '&',"%.2f" % time_dur, '\\\\')
                dontprint = True
            else:
                dontprint = False
            if max_y > best_alt:
                best_alt = max_y
                best_h = i * 0.01 * H
                best_t = where_y
            thistime = time.perf_counter()
            print(' ----------------------------------------------------'+'\n Calculation when p(0)='+str("%.2f"%j)+' and k='+str("%.3f"%(i/100))+' is done.\n Time duration:'+str("%.2f"%(thistime-lasttime))+'\n Time used:'+timeoutput(thistime)+'\n Time last:'+timeoutput(timepredict(thistime-lasttime))+'\n ----------------------------------------------------')
            lasttime = thistime
            #print(best_h)
            h_now[cnt] = max_y
        except:
            thistime = time.perf_counter()
            print(' ----------------------------------------------------'+'\n An error occured when p(0)='+str("%.2f"%j)+' and k='+str("%.3f"%(i/100))+'.\n Time used:'+timeoutput(thistime))
            traceback.print_exc()

    #print(best_t)
    h_best.append(best_h)
    alt_best.append(best_alt)
    t_best.append(best_t * 0.01)
    #print(h_best)
    posi += 1
    #thistime = time.perf_counter()
    #print('----------------------------------'+'\n Calculation when p(0)='+str("%.2f"%j)+' is done.\n Time duration:'+str("%.2f"%(thistime-lasttime))+'\n Time used:'+timeoutput(thistime)+'\n Time last:'+timeoutput(timepredict(thistime-lasttime))+'\n ----------------------------------')
    #lasttime = thistime
    
print(" ----------------------------------------------------")
print(" Caculation DONE!")

hH_best = np.array(h_best) / H
h_eval_b = h_eval / 100

fig3 = plt.figure('fig3',figsize=(4,2), dpi=150)
ax3 = fig3.add_subplot(111)
ax3.plot(h_eval_b,h_now, color='black')
plt.xlabel("初始水高,h(0)/H")
plt.ylabel("最大高度,m")
ax3.xaxis.set_major_locator(MultipleLocator(0.1))
ax3.xaxis.set_minor_locator(MultipleLocator(0.02))
ax3.yaxis.set_major_locator(MultipleLocator(10))
ax3.yaxis.set_minor_locator(MultipleLocator(2))
plt.axis([0,1,0,30])
plt.grid()
plt.tight_layout()
plt.savefig('figs/fig3.png')

plt.show()
