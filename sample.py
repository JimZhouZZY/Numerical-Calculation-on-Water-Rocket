import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
# 物理量赋值
m0,rho,rho0,S0,S,H,p0,g,cd,gamma,k = 0.1,1000,1.2250,4.792e-4,59.447e-4,0.25233,1.01325e5,9.81,0.4,1.4,1/3
ht0,pt0 = H * k, 3*p0
y0 = [ht0, 0, 0]# 初始条件
dydt = [0,0,0]#三个求解的函数，依次分别为：水量、速度、高度
def vabsv(argv):
    return (argv / abs(argv)) if argv != 0 else 1
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
ans = solve_ivp(fun, (0, 10), y0, method='Radau',t_eval = np.arange(0,10,0.01))#求解
plt.plot(ans.t, ans.y[1, :])#绘图
plt.show()
