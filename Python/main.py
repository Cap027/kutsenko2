import math
from CoolProp.CoolProp import PropsSI

class Const:
    g = 9.81
    R = 8.31
    sigma = 5.67e-8
    T00 = 273.15
    M = 29e-3
    p0 = 1013e2

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

        self.tavg = self.t0 + 1#0.5*(self.t0 + self.t1)
        self.G = None
        self.Q = None
        self.alpha = None
        self.t2 = None


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
            if abs(t2 - t2_old) < 0.01:
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
            if abs(Q1 - Q2) < 0.01:
                break
            else:
                tavg += 0.001*(Q2-Q1)

        self.G = G
        self.Q = 0.5*(Q1+Q2)
        self.alpha = alpha
        self.t2 = t2
        self.tavg = tavg


gap = Gap(0, 90)
t1 = [t for t in range(90, 200)]
for t in t1:
    gap.t1 = t
    gap.run()
    print('t1 =', '%.1f' % gap.t1, '\tt2 =', '%.1f' % gap.t2, '\ttavg =', '%.1f' % gap.tavg, '\tQ =', '%.1f' % gap.Q, '\talpha =', '%.1f' % gap.alpha, 'G =', '%.5f' % gap.G)

#tavg = 10
#print(calc_G(tavg))
#print(calc_t2(90, 10))



