import numpy as np
from astroquery.splatalogue import Splatalogue, utils
from astropy import units as u
from scipy import convolve
import pylab as pl
from astropy.table import Column

naclv0 = utils.minimize_table(Splatalogue.query_lines(80*u.GHz, 300*u.GHz, chemical_name='NaCl v = 0', line_lists=['CDMS']))
naclv0.add_column(Column(name='nacl_neighbor4GHz', dtype=int, data=np.zeros(len(naclv0), dtype=int)))
naclv0.add_column(Column(name='nacl_neighbor8_16GHz', dtype=int, data=np.zeros(len(naclv0), dtype=int)))
naclv0.add_column(Column(name='kcl_neighbor4GHz', dtype=int, data=np.zeros(len(naclv0), dtype=int)))
naclv0.add_column(Column(name='kcl_neighbor8_16GHz', dtype=int, data=np.zeros(len(naclv0), dtype=int)))
naclv0.add_column(Column(name='ttl_neighbor', dtype=int, data=np.zeros(len(naclv0), dtype=int)))

kclv0 = utils.minimize_table(Splatalogue.query_lines(80*u.GHz, 300*u.GHz, chemical_name=' KCl v = 0', line_lists=['CDMS']))
kclv0.add_column(Column(name='nacl_neighbor4GHz', dtype=int, data=np.zeros(len(kclv0), dtype=int)))
kclv0.add_column(Column(name='nacl_neighbor8_16GHz', dtype=int, data=np.zeros(len(kclv0), dtype=int)))
kclv0.add_column(Column(name='kcl_neighbor4GHz', dtype=int, data=np.zeros(len(kclv0), dtype=int)))
kclv0.add_column(Column(name='kcl_neighbor8_16GHz', dtype=int, data=np.zeros(len(kclv0), dtype=int)))
kclv0.add_column(Column(name='ttl_neighbor', dtype=int, data=np.zeros(len(kclv0), dtype=int)))

for row in kclv0:
    naclsep = np.abs(row['Freq'] - naclv0['Freq'])
    row['nacl_neighbor4GHz'] = (naclsep < 4).sum()
    row['nacl_neighbor8_16GHz'] = ((naclsep > 8) & (naclsep < 16)).sum()
    kclsep = np.abs(row['Freq'] - kclv0['Freq'])
    row['kcl_neighbor4GHz'] = (kclsep < 4).sum()
    row['kcl_neighbor8_16GHz'] = ((kclsep > 8) & (kclsep < 16)).sum()
    row['ttl_neighbor'] = row['kcl_neighbor8_16GHz'] + row['kcl_neighbor4GHz'] + row['nacl_neighbor8_16GHz'] + row['nacl_neighbor4GHz']


for row in naclv0:
    naclsep = np.abs(row['Freq'] - naclv0['Freq'])
    row['nacl_neighbor4GHz'] = (naclsep < 4).sum()
    row['nacl_neighbor8_16GHz'] = ((naclsep > 8) & (naclsep < 16)).sum()
    kclsep = np.abs(row['Freq'] - kclv0['Freq'])
    row['kcl_neighbor4GHz'] = (kclsep < 4).sum()
    row['kcl_neighbor8_16GHz'] = ((kclsep > 8) & (kclsep < 16)).sum()
    row['ttl_neighbor'] = row['kcl_neighbor8_16GHz'] + row['kcl_neighbor4GHz'] + row['nacl_neighbor8_16GHz'] + row['nacl_neighbor4GHz']


# # clever/fast/not that useful:
# bins = np.arange(80, 300.1, 0.1)
# grid = (bins[1:] + bins[:-1])/2
# naclct,_ = np.histogram(naclv0['Freq'], bins)
# kclct,_ = np.histogram(kclv0['Freq'], bins)
#
# # upper/lower sideband separated by 12 on the inside, 16 outside
# kernel = np.zeros(280)
# kernel[:40] = 1
# kernel[120:160] = 1
# kernel[-40:] = 1
#
# smnaclct = convolve(naclct, kernel, mode='same')
# smkclct = convolve(kclct, kernel, mode='same')
#
#
# #pl.plot(grid, smnaclct)
# #pl.plot(grid, smkclct)
# pl.plot(grid, smkclct*smnaclct)
# pl.plot(grid, kclct)
# pl.plot(grid, naclct)

print('band3')
frqrng = (80, 115)
print(naclv0[(naclv0['Freq'] > frqrng[0]) & (naclv0['Freq'] < frqrng[1])])
print(kclv0[(kclv0['Freq'] > frqrng[0]) & (kclv0['Freq'] < frqrng[1])])

# band 4
print('band4')
frqrng = (125, 165)
print(naclv0[(naclv0['Freq'] > frqrng[0]) & (naclv0['Freq'] < frqrng[1])])
print(kclv0[(kclv0['Freq'] > frqrng[0]) & (kclv0['Freq'] < frqrng[1])])


# band 5
frqrng = (155, 210)
print('band5')
print(naclv0[(naclv0['Freq'] > frqrng[0]) & (naclv0['Freq'] < frqrng[1])])
print(kclv0[(kclv0['Freq'] > frqrng[0]) & (kclv0['Freq'] < frqrng[1])])

print('band6')
frqrng = (210, 275)
print(naclv0[(naclv0['Freq'] > frqrng[0]) & (naclv0['Freq'] < frqrng[1])])
print(kclv0[(kclv0['Freq'] > frqrng[0]) & (kclv0['Freq'] < frqrng[1])])

