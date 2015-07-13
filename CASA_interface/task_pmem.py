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

	l = [nbrpols, FilenameofI, FilenameofQ, FilenameofU, FilenameofV, Filenameofdirtybeam, FilenameofdefaultMap, Filenameofdefaultmask, estimatedspacingflux, conservefluxmode, rmsI, rmsQ, rmsU, rmsV, nbreofiteration, bmaj, bmin, bpa, blcx, blcy, trcx, trcy, accelerationfactor, qfactor, polarisationfactor, outputname, ignoreedgepixels, debug]

  
     

	f = open('mempy_driver.dat','w')	
        for i in range(len(l)):
        	f.write(  str(l[i]) +'\n') 
         
        f.close()
	os.system('/home/quentin/mem_code/pmem-master/pmem')
   
        casalog.origin('pmem')
	
