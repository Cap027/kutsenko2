import math
from CoolProp.CoolProp import PropsSI
import numpy as np
from matplotlib import pyplot as plt

class Const:
    g = 9.81
    R = 8.31
    sigma = 5.67e-8
    T00 = 273.15
    M = 29e-3
    p0 = 1013e2
    t1_max = 90
    t2_max = 60

    @staticmethod
    def tosi(param):
            return param+Const.T00


class Gap:
    def __init__(self, t0, t1):
        self.h = 878e-3
        self.d1 = 585e-3
        self.r1 = self.d1/2
        self.gap = 70e-3
        self.r2 = self.r1 + self.gap
        self.d2 = 2 * self.r2
        self.F = math.pi/4*(self.d2**2 - self.d1**2)
        self.Per = math.pi*(self.d1 + self.d2)
        self.d_h = 4*self.F/self.Per
        self.t0 = t0
        self.t1 = t1
        self.e1 = 0.8
        self.e2 = 0.6
        self.epr = 1 / (1/self.e1 + self.d1/self.d2*(1/self.e2 - 1))

        self.V = 216.5e-3#self.h*np.pi*self.d1**2/4


        self.tavg = self.t0 + 1#0.5*(self.t0 + self.t1)
        self.tout = 2 * self.tavg - self.t0
        self.G = None
        self.Q = None
        self.alpha = None
        self.t2 = None

        self.q = None



    def calc_G(self, nu, tavg):
        return (2*Const.g*self.d_h**(5/4)/0.316/nu**0.25 * (tavg-self.t0)/Const.tosi(self.t0))**(4/7) * Const.p0/Const.R/Const.tosi(tavg) * Const.M*self.F

    #@staticmethod
    def calc_Q_G(self, G, tavg):
        cp = PropsSI('CPMASS', 'T', Const.tosi(tavg), 'P', Const.p0, 'air')
        return 2*G*cp*(tavg-self.t0)

    def calc_Q_T(self, alpha, t2, tavg):
        return math.pi*self.h*self.d1*(alpha*(self.t1-tavg)+self.epr*Const.sigma*(Const.tosi(self.t1)**4 - Const.tosi(t2)**4))

    def calc_alpha(self, Re, tavg):
        lam = PropsSI('CONDUCTIVITY', 'T', Const.tosi(tavg), 'P', Const.p0, 'air')
        Pr = PropsSI('Prandtl', 'T', Const.tosi(tavg), 'P', Const.p0, 'air')
        if Re <= 40:
            C, n, m = 0.76, 0.4, 0.37
        elif 40 < Re <= 1000:
            C, n, m = 0.52, 0.5, 0.37
        elif 1000 < Re <= 2e5:
            C, n, m = 0.26, 0.6, 0.37
        else:
            C, n, m = 0.023, 0.8, 0.4
        Nu = C * Re**n * Pr**m

        return Nu * lam / self.d_h

    def calc_t2(self, alpha, tavg):
        t2_old = tavg
        while True:
            t2 = self.epr*Const.sigma * self.d1 / self.d2 / alpha * (Const.tosi(self.t1)** 4 - Const.tosi(t2_old)** 4) + tavg
            if abs(t2 - t2_old) < 0.001:
                return t2
            t2_old = t2

    def run(self):
        tavg = self.tavg
        while True:
            mu, rho = map(lambda x: PropsSI(x, 'T', Const.tosi(tavg), 'P', Const.p0, 'air'), ('VISCOSITY', 'D'))
            nu = mu/rho
            G = self.calc_G(nu, tavg)
            Q1 = self.calc_Q_G(G, tavg)
            alpha = self.calc_alpha(G*self.d_h/mu/self.F, tavg)
            t2 = self.calc_t2(alpha, tavg)
            Q2 = self.calc_Q_T(alpha, t2, tavg)
            if abs(Q1 - Q2) < 0.001:
                break
            else:
                tavg += 0.001*(Q2-Q1)

        self.G = G
        self.Q = 0.5*(Q1+Q2)
        self.alpha = alpha
        self.t2 = t2
        self.tavg = tavg
        self.tout = 2*tavg - self.t0
        self.q = self.Q / self.V


#gap = Gap(0, 90)
tmin = -20
tmax = 240
temps = np.arange(tmin+1, tmax+1, 0.5)
cases = [Gap(-20, t1) for t1 in temps]
for case in cases:
    case.run()
    print("t1 =", case.t1, "calculated")

t1 = [x.t1 for x in cases]
t2 = [x.t2 for x in cases]
tavg = [x.tavg for x in cases]
tout = [x.tout for x in cases]
Q = [x.Q/1000 for x in cases]

inter_barrel = []
inter_concr = []
inter_barrel.insert(0, min(t1, key=lambda x: abs(x-Const.t1_max)))
inter_barrel.insert(0, Q[t1.index(inter_barrel[0])])
inter_concr.insert(0, min(t2, key=lambda x: abs(x-Const.t2_max)))
inter_concr.insert(0, Q[t2.index(inter_concr[0])])

#print(inter_barrel[0] ,inter_barrel[1])
#print(inter_concr[0], inter_concr[1])
G = [x.G for x in cases]
alpha = [x.alpha for x in cases]
q = [x.q for x in cases]
i = t1.index(inter_barrel[1])
print("Q =", Q[i], "q =", q[i], "t1 =", t1[i], "t2 =", t2[i], "tavg =", tavg[i], "tout =", tout[i], "G =", G[i], "alpha =", alpha[i])
i = t2.index(inter_concr[1])
print("Q =", Q[i], "q =", q[i], "t1 =", t1[i], "t2 =", t2[i], "tavg =", tavg[i], "tout =", tout[i], "G =", G[i], "alpha =", alpha[i])





fig, axs = plt.subplots(nrows=2, ncols=1)
axs[0].plot(Q, t1, label='Бочка', color='black')
axs[0].plot(Q, [Const.t1_max for x in Q], "--", label='Температура вспенивания', color='black')
axs[0].plot(Q, t2, label='Бетон', color='brown')
axs[0].plot(Q, [Const.t2_max for x in Q], "--", label='Температура охрупчивания', color='brown')
axs[1].plot(Q, tavg, color='blue', label='Средняя температура воздуха')
axs[1].plot(Q, tout, color='red', label='Температура воздуха на выходе')

fig.supxlabel("Тепловая мощность, кВт")
fig.supylabel("Температура, град. Цельс")
axs[0].set_xlim(0, 9)
axs[1].set_xlim(0, 9)
axs[0].set_ylim(tmin, 200)
axs[0].yaxis.set_ticks(np.arange(-20, 220+40, 40))
axs[1].set_ylim(tmin, 5)
axs[0].grid()
axs[1].grid()
axs[0].legend()
axs[1].legend()


plt.figure(2)
alpha = [x.alpha for x in cases]
plt.plot(Q, alpha)
plt.xlabel("Тепловая мощность, кВт")
plt.ylabel("Коэффициент теплоотдачи, Вт/$м^2$/К")
plt.xlim(0, 9)
plt.ylim(0, 18)
plt.grid()


plt.figure(3)
G = [x.G for x in cases]
plt.plot(Q, G)
plt.xlabel("Тепловая мощность, кВт")
plt.ylabel("Массовый расход, кг/$м^3$")
plt.xlim(0, 9)
plt.ylim(0, 0.5)
plt.grid()
plt.show()


#t1 = [t for t in range(90, 200)]
#for t in t1:
#    gap.t1 = t
#    gap.run()
#    print('t1 =', '%.1f' % gap.t1, '\tt2 =', '%.1f' % gap.t2, '\ttavg =', '%.1f' % gap.tavg, '\tQ =', '%.1f' % gap.Q, '\talpha =', '%.1f' % gap.alpha, 'G =', '%.5f' % gap.G)

#tavg = 10
#print(calc_G(tavg))
#print(calc_t2(90, 10))



