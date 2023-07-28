import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.integrate import solve_ivp
import numpy as np
import time, traceback
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']

rge = 10
condi = -1e-1

m0 = 0.1
rho = 1000
rho0 = 1.2250
S0 = 4.792e-4
S = 59.447e-4
H = 0.25233
ht0 = 0
p0 = ATM = 1.01325 * 100000
pt0 = 3 * ATM
g = 9.81
cd = 0.4
gamma = 1.4

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

oof = oof0 = 1e5
a = 0
def fun(t, y):
    h,v,pos_y = y[0],y[1],y[2]
    p = pt0 *(((H-ht0)/(H-h))**gamma)
    m = m0 + rho * S * h
    A = (2*S0)/(m*(1-(S0**2/S**2)))
    dydt[2] = v if (pos_y>0 or v>=0) else(0 if (pos_y<=0 or v<0) else 0)
    dydt[1] = (A*(p+rho*g*h-p0)-g-0.5*cd*rho0*S*v*v*vabsv(v)/m)/(1-A*rho*h) if ((pos_y>0 or v>=0) and h>0) else(-g-0.5*cd*rho0*S*v*v*vabsv(v)/ m if ((pos_y>0 or v>=0) and h<=0) else(0 if (pos_y<=0 and v<0) else 0))
    a = dydt[1]
    dydt[0] = -(S0/S)*np.sqrt((2*(p+rho*(g+a)*h-p0))/(rho*(1-(S0**2/S**2)))) if (h>0) else(0 if (h<=0) else 0)
    return dydt

# 初始条件
h_eval = np.arange(1,55,0.01)
p_eval = np.arange(1,10,0.1)

y0 = [0, 0, 0]
dydt = [0,0,0]
N = len(h_eval) * len(p_eval)
h_best = np.zeros((len(p_eval)),dtype=float) + 1 - 1
alt_best = np.zeros((len(p_eval)),dtype=float) + 1 - 1
t_best = np.zeros((len(p_eval)),dtype=float) + 1 - 1
lasttime = pos =0

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
    for i in h_eval:
        N -= 1
        try:
            ht0 = H*i*0.01
            a=0
            oof = 1e5
            y0 = [ht0, 0, 0]
            dydt = [0,0,0]
            yy = solve_ivp(fun, (0, rge), y0, method='Radau',t_eval = np.arange(0,rge,0.01) )
            t = yy.t
            data = yy.y
            max_y = np.amax(data[2])
            where_y = np.argmax(data[2])
            time_dur = find_root(data[2]) * 0.01
            if max_y > best_alt:
                best_alt = max_y
                best_h = i * 0.01 * H
                best_t = where_y
            thistime = time.perf_counter()
            print(' ----------------------------------------------------'+'\n Calculation when p(0)='+str("%.6f"%j)+' and k='+str("%.6f"%(i/100))+' is done.\n Time duration:'+str("%.2f"%(thistime-lasttime))+'\n Time used:'+timeoutput(thistime)+'\n Time last:'+timeoutput(timepredict(thistime-lasttime))+'\n ----------------------------------------------------')
            lasttime = thistime
        except:
            thistime = time.perf_counter()
            print(' ----------------------------------------------------'+'\n An error occured when p(0)='+str("%.6f"%j)+' and k='+str("%.6f"%(i/100))+'.\n Time used:'+timeoutput(thistime))
            traceback.print_exc()
    h_best[pos] = best_h
    alt_best[pos] = best_alt
    t_best[pos] = best_t * 0.01
    pos += 1
    
print(" ----------------------------------------------------")
print(" Caculation DONE!")

hH_best = np.array(h_best) / H

fig4 = plt.figure('fig1',figsize=(4,2), dpi=150)
ax4 = fig4.add_subplot(111)
ax4.plot(p_eval,alt_best, color='black')
ax4.xaxis.set_major_locator(MultipleLocator(1))
ax4.xaxis.set_minor_locator(MultipleLocator(0.2))
ax4.yaxis.set_major_locator(MultipleLocator(20))
ax4.yaxis.set_minor_locator(MultipleLocator(4))
plt.xlabel("初始压强p(0),ATM")
plt.ylabel("最大高度,m")
plt.axis([0,10,0,75])
plt.tight_layout()
plt.grid()
plt.savefig('figs/fig4.png')

fig5 = plt.figure('fig2',figsize=(4,2), dpi=150)
ax5 = fig5.add_subplot(111)
ax5.plot(p_eval,hH_best, color='black')
ax5.xaxis.set_major_locator(MultipleLocator(1))
ax5.xaxis.set_minor_locator(MultipleLocator(0.2))
ax5.yaxis.set_major_locator(MultipleLocator(0.1))
ax5.yaxis.set_minor_locator(MultipleLocator(0.02))
plt.xlabel("初始压强p(0),ATM")
plt.ylabel("最佳装水量k,1")
plt.axis([0,10,0,0.5])
plt.tight_layout()
plt.grid()
plt.savefig('figs/fig5.png')

fig6 = plt.figure('fig3',figsize=(4,2), dpi=150)
ax6 = fig6.add_subplot(111)
ax6.plot(p_eval,t_best, color='black')
ax6.xaxis.set_major_locator(MultipleLocator(1))
ax6.xaxis.set_minor_locator(MultipleLocator(0.2))
ax6.yaxis.set_major_locator(MultipleLocator(1))
ax6.yaxis.set_minor_locator(MultipleLocator(0.2))
plt.xlabel("初始压强p(0),ATM")
plt.ylabel("时间,s")
plt.axis([0,10,0,3.5])
plt.tight_layout()
plt.grid()
plt.savefig('figs/fig6.png')

plt.show()
