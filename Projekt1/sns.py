from datetime import date
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
    
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
data = [2022, 2, 25, 0, 0, 0]  #####################################################################
week, tow = date2tow(data)
# print(week, tow)

nav = navall[0,:]
# print(nav)

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
    # print(tk)  # print(tk/86400)

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
    return np.array([Xk, Yk, Zk])


satelita = satpos(nav, week, tow)
# print(f'{satelita=}')


# Zajęcia 2
def geo2xyz(fi, lam, h, a=6378137, e2=0.00669437999013):
    """funkcja zamienia współrzędne geodezyjne na kartezjańskie fi, lam podajemy w radianach do wzorów, do fcji w deg"""
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    N = a / np.sqrt(1 - e2 * np.sin(fi) ** 2)
    x = (N + h) * np.cos(fi) * np.cos(lam)
    y = (N + h) * np.cos(fi) * np.sin(lam)
    z = (N * (1 - e2) + h) * np.sin(fi)
    return np.array([x, y, z])

def xyz2neu(fi, lam, A, B):
    # def xyz2neu(A:punktodniesienia/startowy, B:koniecwektora):
    """funkcja zamienia wsp kartezjańskie na topocentryczne neu
    A, B reprezentują punkty, A to początek, B to koniec wektora
    A, B są typu np.array i mają 3 współrzędne: x, y, z
    fi, lam to współrzędne punktu A potrzebne do macierzy obrotu"""
    # x, y, z -> north, east, up
    # fi, lambda, lotniska
    # wektor AB
    fi = np.deg2rad(fi)
    lam = np.deg2rad(lam)
    rotation_matrix = np.array([
        [-1*np.sin(fi)*np.cos(lam), -1*np.sin(lam), np.cos(fi)*np.cos(lam)],
        [-1*np.sin(fi)*np.sin(lam), np.cos(lam), np.cos(fi)*np.sin(lam)],
        [np.cos(fi), 0, np.sin(fi)]
    ])
    vector = B - A
    return rotation_matrix.transpose() @ vector

fi, lam, h = 52, 21, 100
odbiornik = geo2xyz(fi, lam, h)  # miejsce obserwacji, niekoniecznie odbiornik
# print(f'{odbiornik=}')
# print(f'{satelita-odbiornik=}')

# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# ax.scatter(3, 2)
# plt.show()

df = pd.DataFrame(navall, columns=['PRN', 'health', 'eccentricity', 'toa', 'oi', 'rora', 'sqrta', 'raaw', 'aop', 'ma', 'af0', 'af1', 'week'])

xsat = []
ysat = []
zsat = []
xodb = []
yodb = []
zodb = []
nn = []
ee = []
uu = []
az = []
el = []

for sat in navall:
    Xs = satpos(sat, week, tow)
    xsat.append(Xs[0])
    ysat.append(Xs[1])
    zsat.append(Xs[2])
    Xr = odbiornik
    xodb.append(Xr[0])
    yodb.append(Xr[1])
    zodb.append(Xr[2])
    # Xsr = Xs - Xr
    neu = xyz2neu(fi, lam, Xr, Xs)
    n, e, u = neu
    nn.append(n)
    ee.append(e)
    uu.append(u)
    Az = np.arctan2(e,n)  # arctan2 ???  # azymut
    ele = np.arcsin(u/np.sqrt(n**2 + e**2 + u**2))  # elewacja
    azst = np.rad2deg(Az)
    elst = np.rad2deg(ele)
    az.append(azst)
    el.append(elst)

df['xsat'] = xsat
df['ysat'] = ysat
df['zsat'] = zsat
df['xodb'] = xodb
df['yodb'] = yodb
df['zodb'] = zodb
df['n'] = nn
df['e'] = ee
df['u'] = uu
df['az'] = az
df['el'] = el


# Zajęcia 3
df.set_index('PRN', inplace=True)
df['distance'] = df.apply(lambda row : np.linalg.norm(np.array([row['xsat'], row['ysat'], row['zsat']]) - np.array([row['xodb'], row['yodb'], row['zodb']])), axis=1)
df['Ax'] = df.apply(lambda row : -(row['xsat'] - row['xodb'])/row['distance'], axis=1)
df['Ay'] = df.apply(lambda row : -(row['ysat'] - row['yodb'])/row['distance'], axis=1)
df['Az'] = df.apply(lambda row : -(row['zsat'] - row['zodb'])/row['distance'], axis=1)

# print(df.info())
maska = 10
maska_df = df[df['el'] > maska]
print(maska_df)
print(maska_df.index.values)
A = maska_df[['Ax', 'Ay', 'Az']].to_numpy()
wiersze = maska_df.shape[0]
ones = np.array([np.repeat(1, wiersze)])
A = np.concatenate((A, ones.T), axis=1)
# print(A)
Q = np.linalg.inv(A.T @ A)
# print(Q)
qx, qy, qz, qt = np.diagonal(Q)
GDOP = np.sqrt(qx+qy+qz+qt)
PDOP = np.sqrt(qx+qy+qz)
TDOP = np.sqrt(qt)
print(f'{GDOP = }\n{PDOP = }\n{TDOP = }')


fi = np.deg2rad(fi)
lam = np.deg2rad(lam)
R = np.array([
        [-np.sin(fi)*np.cos(lam), -np.sin(lam), np.cos(fi)*np.cos(lam)],
        [-np.sin(fi)*np.sin(lam), np.cos(lam), np.cos(fi)*np.sin(lam)],
        [np.cos(fi), 0, np.sin(fi)]])
Qxyz = Q[:-1,:-1]
Qneu = R.T @ Qxyz @ R
# print(Qneu)

qn, qe, qu = np.diagonal(Qneu)
HDOP = np.sqrt(qn+qe)
VDOP = np.sqrt(qu)
print(f'{HDOP = }\n{VDOP = }')


