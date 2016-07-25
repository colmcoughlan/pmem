# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:22:29 2015

@author: quentin
"""
import os
import pickle
from taskinit import *

pmem_path = 'pmem'	# path to the PMEM executable on your system, e.g. /Users/admin/git/pmem/


def pmem(
              npol,
              imap,
              qmap,
              umap,
              vmap,
              dirty_beam,
              default_map,
              mask,
              flux,
              fluxmode,
              inoise,
              qnoise,
              unoise,
              vnoise,
              niter,
              bmaj,
              bmin,
              bpa,
              blcx,
              blcy,
              trcx,
              trcy,
              afactor,
              qfactor,
              polfactor,
              outname,
              edgepixels,
              solver,
              verbose):
	casalog.filter('INFO')
        casalog.origin('pmem')

	l = [npol, imap, qmap, umap, vmap, dirty_beam, default_map, mask, flux, fluxmode, inoise, qnoise, unoise, vnoise, niter, bmaj, bmin, bpa, blcx, blcy, trcx, trcy, afactor, qfactor, polfactor, outname, edgepixels, solver, verbose]

  
     

	with open('pmem.parset', 'w') as f: # create parset file
		for i in range(len(l)):             # write in this file the parameters of the list l (one by line)
			f.write(  str(l[i]) +'\n') 

	os.system(pmem_path) #close the file and execute the pmem code with the right link !
	casalog.filter('INFO')
	casalog.origin('pmem')
'''
	f = open('pmem_log.txt','r')	
	for line in f:
    		casalog.post( line.strip() )
	f.close()
'''
   
