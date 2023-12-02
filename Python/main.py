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
    rho_0 = PropsSI('D', 'T', t0 + 273.15, 'P', p0, 'air')

def calc_dp_h(tavg):
    rho_avg = PropsSI('D', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    return 9.81 * Gap.h * (Gap.rho_0 - rho_avg)

def calc_dp_tr(tavg, G):
    rho_avg = PropsSI('D', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    nu = PropsSI('VISCOSITY', 'T', tavg+273.15, 'P', Gap.p0, 'air')
    u = G / rho_avg / Gap.F
    Re = u * Gap.d_h / nu
    ksi = 0.316/Re**0.25
    return ksi * Gap.h/Gap.d_h * G**2/2/rho_avg/Gap.F**2

def calc_G(tavg=60):
    G = 0.001
    while abs(calc_dp_tr(tavg, G) - calc_dp_h(tavg)) > 0.0001:
        G += 0.001

    return G

def calc_Q(tavg=60):
    cp = PropsSI('CPMASS', 'T', tavg+273.25, 'P', Gap.p0, 'air')
    return 2*calc_G(tavg)*(tavg - Gap.t0)

print(calc_G())
print(calc_Q())

