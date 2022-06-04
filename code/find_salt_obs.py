import radio_beam
from astropy import units as u
from astropy import table
from astroquery.alma import Alma, utils as almautils
from astroquery.splatalogue import Splatalogue, utils as splatutils

# broken
# result = Alma().query(payload={'spatial_resolution':'<0.1', 'science_keyword':
#                                ['High-mass star formation',
#                                 'Disks around high-mass stars']},
#                       public=False)

keywords = ('High-mass star formation',
            'High-mass star formation, Astrochemistry',
            'High-mass star formation, HII regions',
            'High-mass star formation, Inter-Stellar Medium (ISM)/Molecular clouds',
            'High-mass star formation, Intermediate-mass star formation',
            'High-mass star formation, Low-mass star formation',
            'High-mass star formation, Pre-stellar cores, Infra-Red Dark Clouds (IRDC)',
            'Disks around high-mass stars',
            'Disks around high-mass stars, Exo-planets')


kws = "', '".join(keywords)
result = Alma().query_tap(f"select * from ivoa.obscore WHERE spatial_resolution<=0.1 AND science_keyword in ('{kws}') AND science_observation='T'")
print(result)
result = result.to_table()

freqs = {row['member_ous_uid']:
         almautils.parse_frequency_support(row['frequency_support'])
         for row in result}


beams = radio_beam.Beams(u.Quantity(result['spatial_resolution'], u.arcsec))
linesens = u.Quantity([u.Quantity(row['sensitivity_10kms'],
                                  u.mJy/u.beam).to(u.K,
                                                   beam.jtok_equiv(freqs[row['member_ous_uid']][0][0]))
                       for beam, row in zip(beams, result)])
lsens = table.Column(data=linesens, name='Line sensitivity (K)',)
result.add_column(col=lsens)

#result = result[result['Line sensitivity (K)'] < 10]
result = result[result['velocity_resolution'] < 10*u.km/u.s]

salts = {(row['target_name'], row['proposal_id']):
         [Splatalogue.query_lines(fs[0], fs[1], chemical_name='NaCl|KCl', line_lists=['CDMS'])
          for fs in freqs[row['member_ous_uid']]]
         for row in result}
has_salt = [any(salts[(row['target_name'], row['proposal_id'])])
            for row in result]
hsrow = table.Column(data=has_salt, name='Has salt',)
result.add_column(col=hsrow)

def try_minimize(tbl):
    try:
        return splatutils.minimize_table(tbl)
    except Exception as ex:
        #print(ex)
        return tbl

salts = {key: [try_minimize(val) for val in value if any(val) and len(val) > 0]
         for key, value in salts.items()}

print(salts)

# Beuther 2018: 2015.1.00496.S       G351.77-0.54  NONDETECTION OF SALTS?
#
# candidate: 2015.1.01163.S             M17UC1
#
