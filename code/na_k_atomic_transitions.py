########################################################
# Started Logging At: 2021-07-16 10:03:42
########################################################

########################################################
# # Started Logging At: 2021-07-16 10:03:43
########################################################
from astroquery.nist import Nist
from astropy import units as u
import pylab as pl
import string

na = Nist.query(1000*u.AA, 6*u.um, 'Na I')
K = Nist.query(1000*u.AA, 6*u.um, 'K I')

pl.plot(na['Observed'], [int(x.strip('*()blar').lstrip("*()blar")) if any(ii in x for ii in string.digits)  else np.nan for x in na['Rel.']], 's', label='Na I')
pl.plot(K['Observed'], [int(x.strip('*()blar').lstrip("*()blar")) if any(ii in x for ii in string.digits)  else np.nan for x in K['Rel.']], 'o', label='K I')

pl.legend(loc='best')
pl.xlabel("Wavelength (Angstroms)")
pl.ylabel("Relative Intensity from NIST")
pl.savefig('Na_and_K_atomic_transitions.png')
