import math

from CoolProp.CoolProp import PropsSI

class Gap:
    h = 878e-3
    d1 = 585e-3
    r1 = d1/2
    gap = 70e-3
    r2 = r1 + gap
    d2 = 2 * r2
    F = math.pi/4*(d2**2 - d1**2)
    Per = math.pi*(d1 + d2)
    d_h = 4*F/Per
    p0 = 1013e2
    t0 = 0
    t1 = 90
    sigma = 5.67e-8
    e1 = 0.8
    e2 = 0.6
    epr = 1 / (1/e1 + d1/d2*(1/e2 - 1))
    tavg = 1
    rho_0 = PropsSI('D', 'T', t0 + 273.15, 'P', p0, 'air')
    rho_avg = PropsSI('D', 'T', tavg+273.15, 'P', p0, 'air')
    nu = PropsSI('VISCOSITY', 'T', tavg+273.15, 'P', p0, 'air')

def calc_dp_h(tavg):
    #Gap.rho_avg = PropsSI('D', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    return 9.81 * Gap.h * (Gap.rho_0 - Gap.rho_avg)

def calc_dp_tr(tavg, G):
    #rho_avg = PropsSI('D', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    #nu = PropsSI('VISCOSITY', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    u = G / Gap.rho_avg / Gap.F
    Gap.Re = u * Gap.d_h / Gap.nu
    ksi = 0.316/Gap.Re**0.25
    return ksi * Gap.h/Gap.d_h * G**2/2/Gap.rho_avg/Gap.F**2

def calc_G(tavg):
    G = 0.001
    while abs(calc_dp_h(tavg) - calc_dp_tr(tavg, G)) > 0.001:
        G += 0.0001

    return G

def calc_Q(tavg=60):
    cp = PropsSI('CPMASS', 'T', tavg+273.25, 'P', Gap.p0, 'air')
    return 2*calc_G(tavg)*cp*(tavg - Gap.t0)

def calc_alpha(Re, tavg):
    lam = PropsSI('L', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    Pr = PropsSI('Prandtl', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    if Re <= 40:
        C, n, m = 0.76, 0.4, 0.37
    elif Re > 40 and Re <= 1000:
        C, n, m = 0.52, 0.5, 0.37
    elif Re > 1000 and Re <= 2e5:
        C, n, m = 0.26, 0.6, 0.37
    else:
        C, n, m = 0.023, 0.8, 0.4
    Nu = C * Re**n * Pr**m
    return Nu * lam / Gap.d_h

def calc_Q1(t1, t2, tavg):
    G = calc_G(tavg)
    #rho = PropsSI('D', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    nu = PropsSI('VISCOSITY', 'T', tavg + 273.15, 'P', Gap.p0, 'air')
    u = G / Gap.rho_avg / Gap.F
    Re = u * Gap.d_h / nu
    alpha = calc_alpha(Re, tavg)

    return alpha*math.pi*Gap.d1*Gap.h*(t1 - tavg) + Gap.epr*Gap.sigma*math.pi*Gap.h*Gap.d1*((t1+273.15)**4 - (t2+273.15)**4)


def calc_t2(t1, tavg):
    G = calc_G(tavg)
    rho = PropsSI('D', 'T', tavg + 273.15, 'P', Gap.p0, 'air')
    nu = PropsSI('VISCOSITY', 'T', tavg + 273.15, 'P', Gap.p0, 'air')
    u = G / rho / Gap.F
    Re = u * Gap.d_h / nu
    alpha = calc_alpha(Re, tavg)
    t2_old = tavg
    t2 = Gap.epr*Gap.sigma*Gap.d1/Gap.d2/alpha*((t1+273.15)**4 - (t2_old+273.15)**4) + tavg
    while abs(t2-t2_old) > 0.001:
        t2_old = t2
        t2 = Gap.epr*Gap.sigma*Gap.d1/Gap.d2/alpha*((t1+273.15)**4 - (t2+273.15)**4) + tavg

    return t2


tavg = Gap.t0+1
t2 = calc_t2(Gap.t1, tavg)
Q = calc_Q(tavg)
Q1 = calc_Q1(Gap.t1, t2, tavg)

while abs(Q - Q1) > 0.001:
    tavg += 0.1
    t2 = calc_t2(Gap.t1, tavg)
    Q = calc_Q(tavg)
    Q1 = calc_Q1(Gap.t1, t2, tavg)
    print(t2, Q, Q1)

#tavg = 10
#print(calc_G(tavg))
#print(calc_t2(90, 10))

