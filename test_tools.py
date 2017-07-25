import apogee_tools as ap
import os

if __name__ == '__main__':

#--------------------------------------------------------------------
# Downloading

	# ap.download('2M15141711+0044474', type='aspcap')
	# data = ap.Spectrum(id='2M15141711+0044474', type='aspcap')
	# data.plot(items=['spectrum', 'apModel'], save=True)

#--------------------------------------------------------------------
# Searching

	# params = ['TEFF', 'LOGG', 'M_H']
	# ranges = [[-10000,4000], [0,5], [-2,2]]
	# source_table = ap.multiParamSearch(par=params, select=ranges)

#--------------------------------------------------------------------
# Retrieve models

	# ap.download('2M00173608+6642386', type='apstar')
	data = ap.Spectrum(id='2M03290406+3117075', type='aspcap')
	data.plot(items=['spec', 'lines'], yrange=[.6,1.2], save=True, output='/Users/admin/Desktop/2M03290406+3117075.pdf')

	# mdl = ap.getModel(params=[4000,5,0], grid='PHOENIX', xrange=[15200,15800])
	
	# mdl = ap.getModel(params=[3200, 5.0, 0.0], grid='BTSETTLb', xrange=[15200,16940])
	# mdl.plot()
	# chi, spec, mdl = ap.compareSpectra(data, mdl)

	# params = ['TEFF', 'LOGG', 'M_H']
	# ranges = [[-10000,4000], [0,5], [-2,2]]
	# source_table = ap.multiParamSearch(par=params, select=ranges)
