# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:22:29 2015

@author: quentin
"""
import os
import pickle
from taskinit import *


def pmem(
              nbrpols,
              FilenameofI,
              FilenameofQ,
              FilenameofU,
              FilenameofV,
              Filenameofdirtybeam,
              FilenameofdefaultMap,
              Filenameofdefaultmask,
              estimatedspacingflux,
              conservefluxmode,
              rmsI,
              rmsQ,
              rmsU,
              rmsV,
              nbreofiteration,
              bmaj,
              bmin,
              bpa,
              blcx,
              blcy,
              trcx,
              trcy,
              accelerationfactor,
              qfactor,
              polarisationfactor,
              outputname,
              ignoreedgepixels,
              debug):
	casalog.filter('INFO')
        casalog.origin('pmem')

	l = [nbrpols, FilenameofI, FilenameofQ, FilenameofU, FilenameofV, Filenameofdirtybeam, FilenameofdefaultMap, Filenameofdefaultmask, estimatedspacingflux, conservefluxmode, rmsI, rmsQ, rmsU, rmsV, nbreofiteration, bmaj, bmin, bpa, blcx, blcy, trcx, trcy, accelerationfactor, qfactor, polarisationfactor, outputname, ignoreedgepixels, debug]

  
     

	f = open('mempy_driver.dat','w') # create a file named mempy_driver.dat
        for i in range(len(l)):             # write in this file the parameters of the list l (one by line)
        	f.write(  str(l[i]) +'\n') 
         
        f.close()
	os.system('/home/quentin/mem_code/pmem-master/pmem > pmem_log.txt') #close the file and execute the pmem code with the right link !
	casalog.filter('INFO')
	casalog.origin('pmem')
	
	g = open('pmem_log.txt','r')	
	for ligne in g:
    		para = ligne.strip()
    		casalog.post(para)
	g.close()

   
