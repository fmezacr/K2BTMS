"""
Created on Mon Oct  3 11:39:09 2016
@author: fmeza
"""

# Libraries

from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.pyplot as pl
import numpy as np
from astropy.io import ascii
import os
import bisect

from ivs.timeseries import pergrams
from ivs.timeseries import freqanalyse as fa

import matplotlib
matplotlib.rcParams.update({'font.size': 16})

if os.path.exists("output"):
	os.system("rmdir output")

cat0 = 0
cat1 = 0
cat2 = 0
cat3 = 0
cat4 = 0
total = 0 

for filename in os.listdir("input"):
	print "Processing candidate star #{0}...".format(total+1)
	star = filename[0:9]	
	#star = args.star
	os.system("mkdir output/{0}".format(star))	
	
	# Reading Files
	
	hdulist = fits.open("input/{0}".format(filename));
	data = hdulist[1].data
	plotx = 'TIME'
	ploty = 'FLUX'
	
	cols = hdulist[1].columns
	
	xdata = data.field(plotx)
	ydatanan = data.field(ploty)
	
	ydata = np.nan_to_num(ydatanan)	
	
	# ASCII Creation, from .fits to .dat for Light Curve (lc)
	
	ascii.write([xdata, ydata], 'output/{0}/{0}_lc.dat'.format(star),format='no_header')
	
	# Display Light Curve
	
	plt.figure(figsize=(18,6))
	plt.xlabel('Time (d)')
	plt.ylabel('Relative Flux (ppm)')
	plt.title('LIGHT CURVE KASOC - EPIC {0}'.format(star))
	plt.plot(xdata,ydata,'k')
	#plt.xlim(xmin=0, xmax=80)
	plt.savefig('output/{0}/{0}_lc.png'.format(star))
	
	# Calculate Periodogram FOURIER
	
	times,signal = data['TIME'],data['FLUX']
	#times,signal = data['TIME'],data['FLUX']-np.average(data['FLUX'])
	newsignal = np.nan_to_num(signal)
	
	#freqs,ampls = pergrams.scargle(times,signal,f0=0, fn=1, df=0.1, norm='amplitude')
	freqs,ampls = pergrams.scargle(times,newsignal)
	
	# ASCII Creation, from .fits to .dat for Lomb-Scargle periodogram (lsp)
	
	ascii.write([freqs, ampls], 'output/{0}/{0}_lsp.dat'.format(star),format='no_header')
	
	# Display Periodogram FOURIER 
	
	#plt.subplot(212)
	plt.figure(figsize=(18,6))
	plt.xlabel('Frequency [d$^{-1}$]')
	plt.ylabel('Amp')
	plt.title('LOMB SCARGLE PERGRAM - EPIC {0}'.format(star))
	#plt.axis([-10, 100, -0.000002, 0.000002])
	plt.plot(freqs,ampls,'k')
	plt.savefig('output/{0}/{0}_lsp.png'.format(star))

	# Iterative Pre-whitening, SNR > 4
	
	filename = 'output/{0}/{0}_lc.dat'.format(star)
	maxiter = 1000
	#f0=0.
	#fn=20.
	
	data = np.loadtxt(filename)
	times,signal = data[:,0],data[:,1]-np.average(data[:,1])
	#params,model = fa.iterative_prewhitening(times,signal,f0=f0,fn=fn,maxiter=maxiter,threads='max',full_output=True,stopcrit=(fa.stopcrit_scargle_snr,4.,1.))
	#params = fa.iterative_prewhitening(freqs,ampls,f0=2,fn=7,maxiter=6,stopcrit=(fa.stopcrit_scargle_snr,4.,1,))
	params,model = fa.iterative_prewhitening(times,signal,maxiter=maxiter,threads='max',full_output=True,stopcrit=(fa.stopcrit_scargle_snr,4.,1.),prewhiteningorder_snr=True,prewhiteningorder_snr_window=1.)
	#print pl.mlab.rec2txt(params,precision=6)
	
	# Creation of the mainfreqs.txt file.
	
	txt = open("output/{0}/{0}_mainfreqs.txt".format(star),"w")
	
	print >> txt, pl.mlab.rec2txt(params,precision=6)
	
	txt.close()
	
	# Display main frequancies, SNR > 4
	
	plt.figure()
	plt.figure(figsize=(18,6))
	plt.vlines(params['freq'],0,params['ampl'],color='k',linestyle='-')
	plt.xlabel('Frequency [c d$^{-1}$]')
	plt.ylabel('Amplitude [magnitude]')
	plt.title('ALL MODES (SNR > 4) - EPIC {0}'.format(star))
	#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	plt.xlim(0,25)
	#plt.ylim(0,1500)
	plt.savefig('output/{0}/{0}_modes.png'.format(star))

	plt.figure()
	plt.vlines(params['freq'],0,params['ampl'],color='k',linestyle='-')
	plt.xlabel('Frequency [c d$^{-1}$]')
	plt.ylabel('Amplitude [magnitude]')
	plt.title('HO G-MODES - EPIC {0}'.format(star))
	#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	plt.xlim(0,2.99)
	#plt.ylim(0,1500)
	plt.savefig('output/{0}/{0}_hogmodes.png'.format(star))
	
	plt.figure()
	plt.vlines(params['freq'],0,params['ampl'],color='k',linestyle='-')
	plt.xlabel('Frequency [c d$^{-1}$]')
	plt.ylabel('Amplitude [magnitude]')
	plt.title('LO G/P-MODES - EPIC {0}'.format(star))
	#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	plt.xlim(3,5.99)
	#plt.ylim(0,1500)
	plt.savefig('output/{0}/{0}_logpmodes.png'.format(star))
	
	plt.figure()
	plt.vlines(params['freq'],0,params['ampl'],color='k',linestyle='-')
	plt.xlabel('Frequency [c d$^{-1}$]')
	plt.ylabel('Amplitude [magnitude]')
	plt.title('HO P-MODES - EPIC {0}'.format(star))
	#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
	#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
	plt.xlim(6,25)
	#plt.ylim(0,25)
	plt.savefig('output/{0}/{0}_hopmodes.png'.format(star))
	
	plt.close("all")
	
	# Modes counting

	modes = params.freq  

	source = sorted(modes)

	poles = (0,3,6,25)

	output = [source[bisect.bisect_left(source,poles[i]):bisect.bisect_right(source,poles[i+1])] for i in range(len(poles)-1)]

	hog = len(output[0])
	logp = len(output[1])
	hop = len(output[2])
	
	B1 = "0" 
	B2 = "0"
	B3 = "0"
	
	for f in modes:
		if 0 <= f <3 : 
			B1 = "1"
		if 3 <= f < 6:
			B2 = "1"
		if 6 <= f < 25 :
			B3 = "1"
	
	comb = B1 + B2 + B3
	
	if comb == "000" :
		cat = "no-modes"
		cat4 += 1
	elif comb == "100" :
		cat = "SPB"
		cat1 += 1
	elif comb == "011" :
		cat = "beta Cep"
		cat2 += 1
	elif comb == "111" :
		cat = "hybrid"
		cat3 += 1
	else: 
		cat = "OTHERS" 	
		cat0 += 1
		
	os.system("touch output/{0}/{1}".format(star,cat))
	total += 1
	
	# Journal PDF Creation
 	
	from reportlab.pdfgen import canvas
	c = canvas.Canvas('output/{0}/{0}.pdf'.format(star))
	c.drawString(250,790, "EPIC {0}".format(star))
	c.drawImage('output/{0}/{0}_lc.png'.format(star),0,620,600,150)
	c.drawImage('output/{0}/{0}_lsp.png'.format(star),0,460,600,150)
	c.drawImage('output/{0}/{0}_modes.png'.format(star),0,300,600,150)
	c.drawImage('output/{0}/{0}_hogmodes.png'.format(star),60,140,170,150)
	c.drawImage('output/{0}/{0}_logpmodes.png'.format(star),230,140,170,150)
	c.drawImage('output/{0}/{0}_hopmodes.png'.format(star),400,140,170,150)
	c.drawString(90,110, "EPIC {0} has ({1}) HO g-modes, ({2}) LO g- and p-modes and ({3}) HO p-modes.".format(star,hog,logp,hop) )
	c.drawString(90,90, "EPIC {0} ===> {1} candidate.".format(star,cat) )
	c.showPage()
	c.save()
	
#print cat0	
#print cat1
#print cat2	
#print cat3	
#print cat4

cat0p = cat0*100/total
cat1p = cat1*100/total
cat2p = cat2*100/total
cat3p = cat3*100/total
cat4p = cat4*100/total

print "---------------------------------------"	
print "RESULTS"	
print "Total of others: {0}%".format(cat0p)	
print "Total of SPB: {0}%".format(cat1p)	
print "Total of beta Cep: {0}%".format(cat2p)	
print "Total of hybrids: {0}%".format(cat3p)	
print "Total of no-modes: {0}%".format(cat4p)	
