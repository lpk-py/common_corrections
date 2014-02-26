import numpy as np
from py2mat_mod import py2mat 

fio = open('./ell_ccor.dataset1_wri_example')
fi = fio.readlines()
data = np.loadtxt('ell_ccor.dataset1_wri_example', delimiter=',')

lat = []; lon = []; rtarg = []
ecorr = []; tau = []; telev = []
for i in range(len(data)):
    lat.append(data[i][0])
    lon.append(data[i][1])
    rtarg.append(data[i][2])
    ecorr.append(data[i][3])
    tau.append(data[i][4])
    telev.append(data[i][5])

py2mat(lat, 'lat', 'lat')
py2mat(lon, 'lon', 'lon')
py2mat(rtarg, 'rtarg', 'rtarg')
py2mat(ecorr, 'ecorr', 'ecorr')
py2mat(tau, 'tau', 'tau')
py2mat(telev, 'telev', 'telev')
