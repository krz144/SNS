import numpy as np
from datetime import date

def readrnxnav(file):
    m=1
    nav=np.zeros((2000,37))
    inav=np.zeros((2000))
    n=-1
    with open(file, "r") as f:
        for s in f:                
            answer = s.find('END OF HEADER') # skip header
            if answer != -1:
                break
        for s in f:
            s = s.replace('D', 'E')
            if m==1:
                prn=int(s2n(s,1,2))
                a = np.empty((1,6))
                a[:] = np.NaN
                a[0,0:6]=np.array(s2e(s,4,23))
            else:
                a=np.append(a,s2n(s,4,19))
            for x in range (3):
                p=23+x*19
                a=np.append(a,s2n(s,p,19))
            if m<8:
                m+=1
            else:
                n+=1
                nav[n,:]=a
                inav[n]=prn
                m=1
        nav=nav[0:n+1,:]
        inav=inav[0:n+1]
        inav = inav.astype(int)
    f.close()
    return nav, inav

def readrnxobs(file, time_start, time_end, GNSS = 'G'):
    with open(file, "r") as f: 
        for s in f:
            label = s[59:]
            if label.find('SYS / # / OBS TYPES') == 1:
                if s[0] == GNSS:
                    p = 7
                    types_header = []
                    for i in range(int(s[4:4+2])):
                        if p > 58:
                            p = 7
                            s = next(f)
                        types_header.append(s[p:p+3])
                        p += 4
                
            elif label.find('END OF HEADER') == 1:
                break
            types_of_obs = ['C1C']
        ind = np.zeros((len(types_header)))
        for n in range(len(types_of_obs)):
            i=(types_header.index(types_of_obs[n])) if types_of_obs[n] in types_header else -1#np.empty((0))
            if i>-1:
                ind[i]=n+1
        
        obs = np.zeros((150000, len(types_of_obs)))*np.nan
        iobs = np.zeros((150000, 3))
        n = 0
        for s in f:
            label = s[0]
            if label == '>':
                epoch = s2e(s,2,29)
                y = epoch[0]
                # tt = (date.toordinal(date(epoch[0],epoch[1],epoch[2]))+366-t0)*86400+np.dot((epoch[3:6]), ([3600,60,1])) + 6*86400
                tt = date2tow(epoch)[1] - date2tow(epoch)[2] * 86400
                if tt > (date2tow(time_end)[1] - date2tow(epoch)[2] * 86400):
                    break
                else:
                    flag = int(round(tt))>=(date2tow(time_start)[1] - date2tow(time_start)[2] * 86400)
                if flag:
                    number_of_all_sats = int(s[33:33+2])
                    iobs[n+np.arange(0,number_of_all_sats),1] = tt
                    iobs[n+np.arange(0,number_of_all_sats),2] = date2tow(epoch)[1]
                    for sat in range(number_of_all_sats):
                        s = next(f)
                        p = 3
                        if s[0] == GNSS:
                            for i in range(len(types_header)):
                                if ind[i] != 0:
                                    obs[n+sat, int(ind[i] - 1)] = s2n(s,p,16)
                                    iobs[n+sat,0] = s2n(s,1,2)
                                p+=16
                    n += number_of_all_sats
        obs = obs[0:n, :]
        iobs = iobs[0:n,:]
        obs = np.delete(obs,iobs[:,0]==0, axis=0)
        iobs = np.delete(iobs,iobs[:,0]==0, axis=0)
        f.close()
        iobs = iobs.astype(int)
    return obs, iobs

def s2e(s,p,n):
    epoch = [int(s[p:p+4]), int(s[p+5:p+5+2]), int(s[p+8:p+8+2]), int(s[p+11:p+11+2]), int(s[p+14:p+14+2]), float(s[p+17:n])]
    return epoch 

def date2tow(data):    
    dday=date.toordinal(date(data[0],data[1],data[2])) - (date.toordinal(date(1980,1,6)))  # ktorys kwietnia 2019 w poprzedniej date2tow, tutaj liczymy od 1980 roku, wczesniej byl rollover
    week = dday//7
    dow = dday%7
    tow = dow * 86400 + data[3] * 3600 + data[4] * 60 + data[5]
    return week, tow, dow

def s2n(s,p,n):
    a = s[p:p+n]
    if (not (a and not a.isspace())):
        a = np.nan
    else:
        a = float(a)        
    return a



plik = 'WROC00POL_R_20220800000_01D_GN.rnx'
nav, inav = readrnxnav(plik)
# print(nav, inav)

sat = 1
ind = inav==sat
nav1 = nav[ind,:]
# print(ind)
# print(nav1[:,:])

t = 86400 + 0.5*3600 # + 12 * 3600  # dla godziny 12 wybieramy ramke


# poczatek funkcji
# arg wejsciowe:
    # tablica nav dla wybranego satelity
    # epoka, na którą chcemy obliczyć współrzędne t
    # tydzień - week (nie musi być to tydzień.. wtf)
    # 

# Omega zero - nav19, itd/ jak w poprzedniej funkcji sat 
# z nav0 17 element to bedzie toe     # t-toe = tk

# week = 2202
# week = 2202

# week = nav1[:,27]  # week jest zapisany gdzies w danych
week = nav1[:,27][0]  # week jest zapisany gdzies w danych
# print(week)

# DATA = [2022,4,7,0,0,0]
# print(date2tow(DATA))


t = t + week * 86400*7
toe_all = nav1[:,17] + nav1[:,27] * 86400*7  # uodpornienie na dane z nd lub soboty wieczor + nav1[:,27] * 86400*7

roznica = t - toe_all
ind_t = np.argmin(abs(roznica))
print(f'{ind_t=}')

nav0 = nav1[ind_t, :]
# print(nav1)
print(nav0)

# WYNIKI na teams DLA SATELITY nr 1 NA GODZINE 00:30
# wyjscie : X,Y,Z, dts  # delta ts cos z zegarem satelity?

def satpos(nav, week, tow):
    mi = 3.986005e14 
    omge = 7.2921151467e-5
    # parametry zegara
    toc = nav[0:5]
    toe = nav[17]
    af0 = nav[6]
    af1 = nav[7]
    af2 = nav[8]
    week = nav[27]
    # elementy orbity keplerowskiej
    sqrta = nav[16]
    e = nav[14]
    i0 = nav[21]
    Omega_zero = nav[19]
    omega = nav[23]
    M0 = nav[12]
    # parametry perturbacyjne
    delta_n = nav[11]
    Omega_dot = nav[24]
    IDOT = nav[25]
    Cuc = nav[13]
    Cus = nav[15]
    Cic = nav[18]
    Cis = nav[20]
    Crc = nav[22]
    Crs = nav[10]

    tk = tow - toe
    # print(f'{tk=}')

    a = sqrta**2
    n0 = np.sqrt(mi/(a**3))
    n = n0 + delta_n
    Mk = M0 + n * tk
    Epop = Mk
    while True:
        E = Mk + e * np.sin(Epop)
        if (abs(E-Epop)<10**(-12)):
            break
        Epop = E
    Ek = E
    vk = np.arctan2(np.sqrt(1 - e**2)*np.sin(Ek) , np.cos(Ek)-e)
    Fik = vk + omega
    duk = Cus * np.sin(2*Fik) + Cuc * np.cos(2*Fik)
    drk = Crs * np.sin(2*Fik) + Crc * np.cos(2*Fik)
    dik = Cis * np.sin(2*Fik) + Cic * np.cos(2*Fik)
    uk = Fik + duk
    rk = a*(1 - e*np.cos(Ek)) + drk
    ik = i0 + IDOT * tk + dik
    xk = rk  * np.cos(uk)
    yk = rk * np.sin(uk)
    Omega_k = Omega_zero + (Omega_dot - omge) * tk - omge * toe
    Xk = xk * np.cos(Omega_k) - yk * np.cos(ik) * np.sin(Omega_k)
    Yk = xk * np.sin(Omega_k) + yk * np.cos(ik) * np.cos(Omega_k)
    Zk = yk * np.sin(ik)
    # print(f'{Xk=} {Yk=} {Zk}')  # Współrzędne satelity w układzie ECEF

    # obliczenie błędu synchronizacji zegara satelity
    delta_t_s = af0 + af1*(t-toe) + af2*(t-toe)**2
    # print(f'{delta_t_s=}')
    # obliczenie poprawki relatywistycznej i dodanie jej do delta_t_s
    c = 299792458.0
    delta_t_rel = ((-2*np.sqrt(mi))/(c**2)) * e*sqrta*np.sin(Ek)
    # print(f'{delta_t_rel=}')
    delta_tsrel = delta_t_s + delta_t_rel
    # print(f'{delta_tsrel=}')

    return(Xk, Yk, Zk, delta_tsrel)


DATA = [2022,3,21,0,30,0]
week, tow, dow = date2tow(DATA)

# print(week, tow, dow)
X, Y, Z, dtrel = satpos(nav0, week, tow)
print(X, Y, Z, dtrel)






