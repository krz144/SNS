import numpy as np


def md(el):
    # el = angle in RADIANS - przyjmuje np.sin
    # el - PODAJEMY W STOPNIACH DO PIERWIASTKA i do funkcji
    # caly pierwiastek na radiany przeliczyc
    return 1/np.sin(np.radians(np.sqrt(el**2 + 6.25)))


def mw(el):
    return 1/np.sin(np.radians(np.sqrt(el**2 + 2.25)))


# wszystkie argumenty: h(GRS80), el_sat
def poprawka_toposferyczna(h_el, el_sat, N=40):
    # Model Saastamoinena
    h = h_el - N

    p0 = 1013.25  # hPA
    t0 = 291.15  # K
    Rh0 = 0.5  # = 50%

    # p-cisnienie t-temp Rh-wilgotnosc wzgledna e-cisnienie pary wodnej
    p = p0 * (1 - 0.0000226*h)**(5.225)
    t = t0 - 0.0065*h
    Rh = Rh0 * np.exp(-0.0006396)
    e = 6.11 * Rh * 10**((7.5*(t-273.15))/(t-35.85))

    deltaT_d0 = 0.002277*p
    deltaT_w0 = 0.002277*(1255/t + 0.05)*e

    deltaT = md(el_sat) * deltaT_d0 + mw(el_sat) * deltaT_w0

    return deltaT


if __name__ == '__main__':
    print(poprawka_toposferyczna(180, el_sat=30, N=40))
