from datetime import date
import numpy as np
    
def read_yuma(almanac_file):
    ''' 
    Reading and parsing YUMA asci format
    INPUT:
        Almanac: YUMA format 
    OUTPUT:
        almanac_data -  type list of list [strings value], number of lists is equal to number of satellite
                        one list contain satellites according to the order:         
                        ['SV ID', 'Health', 'Eccentricity', 'Time of Applicability(s)', 'Inclination(rad)', 
                        'Rate of Right Ascen(r/s)', 'SQRT(A)  (m 1/2)', 'Right Ascen at Week(rad)', 
                        'Argument of Perigee(rad)', 'Mean Anom(rad)', 'Af0(s)', 'Af1(s/s)', 'Week no']
        
    '''
    
    if almanac_file:
        alm = open(almanac_file)
        
        alm_lines = alm.readlines()
        all_sat = []
        for idx, value in enumerate(alm_lines):
            # print(idx, value)
            
            if value[0:3]=='ID:':
                one_sat_block = alm_lines[idx:idx+13]
                one_sat = []
                for line in one_sat_block:
                    data = line.split(':')
                    one_sat.append(float(data[1].strip()))
                all_sat.append(one_sat)
        alm.close()
        all_sat = np.array(all_sat)
        return (all_sat)  


def date2tow(data):
    """
    Parameters
    data : data -- list [year,month,day,hour,minute,second]
    Returns
    week : GPS week, for the second rollover, in range 0-1023
    tow : second of week.
    """
    # difference of days
    dd = date.toordinal(date(data[0], data[1], data[2])) - date.toordinal(date(2019, 4, 7))    
    # week number
    week = dd // 7
    #day of week
    dow = dd % 7
    # time of week
    tow = dow * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
    return week, tow


navall = read_yuma('yumagood.txt')  # 'yumagood.txt'
data = [2022, 2, 25, 0, 0, 0]
week, tow = date2tow(data)
print(week, tow)

nav = navall[0,:]
print(nav)

def satpos(nav, week, tow):
    """jedna duża funkcja, algorytm na TEAMS, zwracać ma x, y, z"""
    # stałe
    mi = 3.986005e14  # prędkość kątowa obrotu Ziemii - to, czy omge?
    omge = 7.2921151467e-5

    prn = nav[0]
    e = nav[2]
    toa = nav[3]
    i = nav[4]
    Omega_dot = nav[5]
    sqrta = nav[6]
    Omega = nav[7]
    omega = nav[8]
    M0 = nav[9]
    gps_week = nav[12] 

    t = week * 7 * 86400 + tow
    toa_weeks = gps_week * 7 * 86400 + toa

    # tk = tow - toa
    tk = t - toa_weeks
    print(tk)  # print(tk/86400)

    # Krok 2 algorytmu
    a = sqrta**2
    n = np.sqrt(mi/(a**3))
    Mk = M0 + n * tk
    # Ek = Mk + e * np.sin(EK)  # rozwiazanie iteracyjne
    Epop = Mk
    while True:
        E = Mk + e * np.sin(Epop)
        # print(E)
        if (abs(E-Epop)<10**(-12)):
            break
        Epop = E
    Ek = E
    # atan2 z math lub arctan2 z numpy
    vk = np.arctan2(np.sqrt(1 - e**2)*np.sin(Ek) , np.cos(Ek)-e)
    Fik = vk + omega
    rk = a * (1 - e * np.cos(Ek))
    xk = rk * np.cos(Fik)
    yk = rk * np.sin(Fik)
    # na dole toa czy toa_weeks?
    Omega_k = Omega + (Omega_dot - omge) * tk - omge * toa  # skąd omega_E?
    Xk = xk * np.cos(Omega_k) - yk * np.cos(i) * np.sin(Omega_k)
    Yk = xk * np.sin(Omega_k) + yk * np.cos(i) * np.cos(Omega_k)
    Zk = yk * np.sin(i)
    return (Xk, Yk, Zk)


wynik = satpos(nav, week, tow)
print(wynik)




