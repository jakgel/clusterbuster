#!usr/Workspace/Promotion0eps/Praktikum
# load each file ---> get images in NRSS, VLSS and FIRST images


#====  InitialiseClustersSQL_v1.1.py =====#
# Writes an .db file for galaxy clusters out of an .cv/dat? with given structure
# Authors: Jakob Gelszinnis 
# v1.0 (Mid 2014)   - Basic utilities
# v1.1 (Jan 2015)   - Fixed many issues
# v1.2 (March 2015) - directly out of the cluster csv, possibility as a subfunction of ExtractRelic.py

# BUG flexfile is not working 
#usage
# e.g. python GetFITSImages_v1.1.py  -i SortedClusters.dat -fs 1


import urllib, urllib2
import os
import ephem
from cosmocalc import cosmocalc
import argparse
import csv
from pyutil.ClusterBuster_pythonClass import *
from astropy.coordinates import SkyCoord
from astropy import units as u
  
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end].replace('"','')
    except ValueError:
        return ""

# This kind of conversion seems much to bloated     
def decdeg2hms(dd):
   is_positive = dd >= 0
   dd = abs(dd/15)
   minutes,seconds = divmod(dd*3600,60)   
   degrees,minutes = divmod(minutes,60)
   degrees = degrees if is_positive else -degrees
   return (str(int(degrees)) + ':' + str(int(minutes)) + ':' + str(seconds))


    
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", type=str, help="file with coordinates", default='ClusterList/Cluster.csv')
parser.add_argument("-fs", "--flexscale", type=int, help="apply a flexible scale for the images downloaded?", default=0)
args = parser.parse_args()


#print >> file, "   #detid     RA           DEC       z       L     Scale   One MParsec in Arcsec"
infile = args.infile

          
ClList = []

def Get_fits(ClusterFile):
  
  
    flexscale   = True #Because of the BUG!
    surveys = ['TGSS_alt'] #NVSS,WENNS,FIRST,TGSS
  

    print '###==== Step 0: Load pointings for TGSS  ====###'
    pointings = []
    for line in open( 'Surveys/TGSS_pointings.dat' ) :

        token = line.rstrip('\n').split()
        if len(token)<2 or token[0][0]=='#' :
          continue
        pointings.append ( (token[0], SkyCoord(ra=float(token[1])*u.degree, dec=float(token[2])*u.degree) ) )
    print 'Number of TGSS pointings:', len(pointings)
  
    print '###==== Step 1: Load cluster file and get parameters  ====###'

    print 'Opening the file "' + infile + '" to get the clusters for search in image databases'
    with open(ClusterFile, 'rb') as csvfile:
      spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
      for CL in spamreader:
	if CL[0] and not CL[0][0] == '#' and  CL[0] not in [o.name for o in ClList]: 
	
	  Cl_name = CL[0]
	  if CL[12] in ['TRUE'] and Cl_name not in ['']:   #If Flag_INCLUDE==TRUE   ; and Cl_name not in ['1RXS_J060313.4+421231']
	    
	    RA     = float(CL[2]) # d.XXXX
	    Dec    = float(CL[3]) # d.XXXX
	  
	    try:
	      z = float(CL[4])
	    except :
	      z = 0.
	      
	    cosmoPara   =  cosmocalc(z, H0=70, WM=0.29, WV=0.71)
	    D_lum       =  cosmoPara['DL_cm']                     # cm 
	    plsc        =  cosmoPara['PS_kpc']                    # kpc/''
        
	    #create Class object
	    GCl     = galaxycluster(name=Cl_name, RA=RA, Dec=Dec, z=z)
	    ClList.append(GCl)
        
	    print '###==== Step 2: Fork fits files from servers for %s ====###'   % (Cl_name)
	    

	    det_id  = Cl_name
	    RA_h    = str(RA/15) # h.XXXX
	    
	    RA_hms  = str( ephem.hours(   str(decdeg2hms(float(RA)) ) ) ).replace( ':', ' ') # hh mm ss
	    Dec_dms = str( ephem.degrees( str(float(Dec)) ) ).replace( ':', ' ')             # dd mm ss
	    
	    RA_hms_p  = str( ephem.hours(   str(decdeg2hms(float(RA)) ) ) ).replace( ':', '+') # hh+mm+ss
	    Dec_dms_p = str( ephem.degrees( str(float(Dec)) ) ).replace( ':', '+')             # dd+mm+ss
	    
	    onempc  = 1e3/plsc       #Imagesize = 30  # in Bogensekunden
	
	    if onempc > 700 and flexscale:
		imN    = ['2+2','15+15']
		imV    = ['3+3','22.5+22.5']
	    else:
		imN    = ['1+1','7.5+7.5']
		imV    = ['2+2','15+15']
		
	    adressnvss  = 'http://www.cv.nrao.edu/cgi-bin/postage.pl?ObjName=NVSS-'     + det_id + '&RA=' + RA_hms_p + '&Dec=' + Dec_dms_p + '&Size='+imN[0]+'&Cells='+imN[1]+'&Equinox=2000&Type=application/octet-stream'
	    adressfirst = 'http://third.ucllnl.org/cgi-bin/firstcutout?RA=' + RA_h + '&Dec=' + str(Dec) + '&ImageSize=30&MaxInt=1&ImageType=FITSImage' #Image Size= 1 equals 1 arcmin
	    adressvlss  = 'http://www.cv.nrao.edu/cgi-bin/newVLSSpostage.pl?ObjName=VLSS' + det_id + '&RA=' + RA_hms_p + '&Dec=' + Dec_dms_p + '&Size='+imV[0]+'&Cells='+imV[1]+'&Equinox=2000&Type=application/octet-stream'   #http://www.cv.nrao.edu/vlss/VLSSpostage.shtml
	    
	    #=== NVSS
	    if 'NVSS' in surveys:
	      print adressnvss
	      for h in range (0,5):
		try:
		  web_file = urllib2.urlopen(adressnvss, timeout = 18)
		  out_file = open('ClusterRegionImages/NVSS/NVSS-' + det_id + '.fits', 'w')
		  out_file.write(web_file.read())
		  out_file.close()
		except:
		  continue
		else: 
		  break        

	    #=== FIRST
	    if 'FIRST' in surveys:
	      print adressfirst 
	      for h in range (0,5):
		try:    
		  web_file = urllib2.urlopen(adressfirst, timeout = 18)
		  out_file = open('ClusterRegionImages/FIRST/FIRST-' + det_id + '.fits', 'w')
		  out_file.write(web_file.read())
		  out_file.close()
		except:
		  continue
		else: 
		  break
	    
	    
	    #=== VLSS
	    if 'VLSS' in surveys:
	      print adressvlss
	      for h in range (0,5):
		try:    
		  web_file = urllib2.urlopen(adressvlss, timeout = 18)
		  out_file = open('ClusterRegionImages/VLSSr/VLSSr-' + det_id + '.fits', 'w')
		  out_file.write(web_file.read())
		  out_file.close()
		except:
		  continue
		else: 
		  break

	      
	    #=== WENSS    
	    if 'WENNS' in surveys:
	      #Webside is easily callable, mayby only one request neccesairy
	      url = 'http://www.astron.nl/wow'
	      values = {'POS'     : RA_hms + ' ' + Dec_dms,
			'SIZE'    : '2.5, 2.5',
			'Equinox' : 'J2000',
			'cutout'  :  1,
			'surv_id' :  1}
			

	      for h in range (0,1):   #range (0,5)
		try:    
		
		    data = urllib.urlencode(values)
		    #print 'Data:',data 
		    queurl = url + '/testcode.php?' + data  
		    print 'WENSSurl:', queurl
		    response  = urllib2.urlopen(queurl)
		    print 'Response:',response    
		    
		    #print response.read()
		    #print 'Response.geturl():',response.geturl()
		    #print response.info()
		    
		    regionurl = ''
		    for line in response:
		      #print line.rstrip()
		      
		      #== Not used until a Region larger than 1*1 deg can get choosen
		      #if not line.find('.fits.gz">cut') == -1:
			    ##print line
			    #substr = find_between( line , 'class="testo3"><a ', '</a></span><br/>') 
			    ##print  'Found in between', substr
			    #regionurl = url + '/' + find_between( substr, 'href=', '.fits.gz>') + ".fits.gz"
			    #print  'Regionurl:',regionurl
		      if not line.find('.GZ</a>') == -1:
			    substr = find_between( line , 'a href="survey1', '</a></span><br/>') 
			    #print 'Substr', substr
			    regionurl = url + '/survey1/' + find_between( substr, '/', '.GZ>') + ".GZ"
			    print  'WENSSRegionurl:',regionurl
	      
	      
		    if not regionurl == '':
		      web_file = urllib2.urlopen(regionurl, timeout = 18)
		      filename = 'ClusterRegionImages/WENSS/WENSS-' + det_id + '.fits.gz'
		      out_file = open(filename, 'w')
		      out_file.write(web_file.read())
		      out_file.close()
		      os.system('gzip -d ' + filename) 
		    
		      print 'WENSS ok for', det_id, '\n'       
		    
		    
		except:
		    continue
		else: 
		    break
		  
		
	    #=== TGSS 
	    if 'TGSS' in surveys:
		url = 'http://tgss.ncra.tifr.res.in/cgi-bin/tgss_postage.cgi'
		values = {'raval' : RA_hms,
			  'decval': Dec_dms,
			  'szval' : '1.0',
			  'szunit': 'deg',
			  'fmtval': 'fits'}
		    
		for h in range (0,1):   #range (0,5)
		  try:    
		      data = urllib.urlencode(values)
		      print 'Data:',data    
		      req = urllib2.Request(url, data)
		      print 'Request:',req  

		      response  = urllib2.urlopen(req)
		      ##print response.read()
		      print 'Response.geturl():',response.geturl()
		      ##print response.info()
		      
		      regionurl = ''
		      for line in response:
			#print line.rstrip()
			if not line.find('<a href="opendat/tgss_') == -1:
			      regionurl = 'http://tgss.ncra.tifr.res.in/' + find_between( line , "<a href=", ">")
			      print  'TGSSurl:',regionurl
		
		      if not regionurl == '':
			web_file = urllib2.urlopen(regionurl, timeout = 18)
			out_file = open('ClusterRegionImages/TGSS/TGSS-' + det_id + '.fits', 'w')
			out_file.write(web_file.read())
			out_file.close()

			print 'TGSS ok for', det_id, '\n'

		      
		  except:
		      continue
		  else: 
		      break
		    
	    #=== TGSS altenrative
	    if 'TGSS_alt' in surveys:
	        
	        try:
		  size = os.path.getsize(os.path.join(os.path.abspath("."), 'Images_TGSS/TGSS-' + det_id + '.fits') )
		except: 
		  size = -1
	        
		if size < 10000:
		    if size > 0:
		       os.remove(os.path.join(os.path.abspath("."), 'Images_TGSS/TGSS-' + det_id + '.fits'))
		    print 'Target ' + det_id + ' was not found in the TGSS-catalogue. It is going to be downloaded.'
	        
		    ClusterPos = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree)
		    dist = [p[1].separation(ClusterPos)  for p in pointings]
		    mindist =  min(dist)
		    print 'TGSS-field minimal distance:', mindist
		    index   = dist.index(mindist)
		    cutname = pointings[index][0]
		    
		    url = 'http://vo.astron.nl/tgssadr/q_fits/imgs/form' # http://vo.astron.nl/tgssadr/q_fits/imgs/form/login?nextURL=http%3A%2F%2Fvo.astron.nl%2Ftgssadr%2Fq_fits%2Fimgs%2Fform
		    values = {'hPOS' : '%8.4f , %8.4f' % (RA, Dec),
			      'hSIZE': 0.5,
			      '_FORMAT':"HTML",
			      'submit':'submit'
			      }
			#_queryForm"
			
		    
			
		    for h in range (0,1):   #range (0,5)
		      try:    
			  #data = urllib.urlencode(values)
			  #print 'Data:',data    
			  #req = urllib2.Request(url, data)
			  #print 'Request:',req  

			  #response  = urllib2.urlopen(req)
			  ###print response.read()
			  #print 'Response.geturl():',response.geturl()
			  ##print response.info()
			  
			  #regionurl = ''
			  #for line in response:
			    ##print line.rstrip()
			    #if not line.find('<a href="opendat/tgss_') == -1:
				  #regionurl = 'http://tgss.ncra.tifr.res.in/' + find_between( line , "<a href=", ">")
				  #print  'TGSSurl:',regionurl
		    
			  regionurl = 'http://vo.astron.nl/getproduct/tgssadr/fits/TGSSADR_%s_5x5.MOSAIC.FITS' % (cutname) #<a href="http://vo.astron.nl/getproduct/tgssadr/fits/TGSSADR_R43D11_5x5.MOSAIC.FITS" class="productlink">TGSSADR_R43D11_5x5.MOSAIC.FITS</a>    
			  print regionurl
			  if not regionurl == '':
			    web_file = urllib2.urlopen(regionurl, timeout = 18)
			    out_file = open('Images_TGSS/TGSS-' + det_id + '.fits', 'w')
			    #out_file = open('ClusterRegionImages/TGSS/TGSS-' + det_id + '.fits', 'w')
			    out_file.write(web_file.read())
			    out_file.close()

			    print 'TGSS ok for', det_id, '\n'

			  
		      except:
			  continue
		      else: 
			  break
		else:
		   print 'Target ' + det_id + ' already in the TGSS-catalogue'
		  
	  #import re

	  #s = 'asdf=5;iwantthis123jasd'
	  #result = re.search('asdf=5;(.*)123jasd', s)
	  #print result.group(1)


	  #import re
	  #r = re.compile('Master(.*?)thon')
	  #m = r.search(str1)
	  #if m:
	      #lyrics = m.group(1)
	      
	      
	  #from urllib2 import Request, urlopen, URLError
	  #req = Request(someurl)
	  #try:
	      #response = urlopen(req)
	  #except URLError as e:
	      #if hasattr(e, 'reason'):
		  #print 'We failed to reach a server.'
		  #print 'Reason: ', e.reason
	      #elif hasattr(e, 'code'):
		  #print 'The server couldn\'t fulfill the request.'
		  #print 'Error code: ', e.code
	  #else:
	      ## everything is fine
		  

		  
	    #Somehow Bugy
	    #try:
	      #if os.path.getsize(os.path.join(os.path.abspath("."), 'ClusterRegionImages/NVSS/NVSS-' + det_id + '.fits') ) < 500:
		  #os.remove(os.path.join(os.path.abspath("."), 'ClusterRegionImages/NVSS/NVSS-' + det_id + '.fits'))
		  #print 'Target ' + det_id + ' was not found in the NVSS-catalogue'  
	    #except:
		  #print 'Target ' + det_id + ' was not found in the NVSS-catalogue, there was also no junk file'  
	    #else:
		  #print 'Target ' + det_id + ' was not found in the NVSS-catalogue, there was also no junk file!'  
		    
	    #try:
	      #if os.path.getsize(os.path.join(os.path.abspath("."), 'ClusterRegionImages/FIRST/FIRST-' + det_id + '.fits') ) < 10000:
		#os.remove(os.path.join(os.path.abspath("."), 'ClusterRegionImages/FIRST/FIRST-' + det_id + '.fits'))
		#print 'Target ' + det_id + ' was not found in the FIRST-catalogue'
	    #except:
		  #print 'Target ' + det_id + ' was not found in the FIRST-catalogue, there was also no junk file'  
	    #else:
		  #print 'Target ' + det_id + ' was not found in the FIRST-catalogue, there was also no junk file'  
	      
	    #try:
	      #if os.path.getsize(os.path.join(os.path.abspath("."), 'ClusterRegionImages/VLSSr/VLSSr-' + det_id + '.fits') ) < 500:
		#os.remove(os.path.join(os.path.abspath("."), 'ClusterRegionImages/VLSSr/VLSSr-' + det_id + '.fits'))
		#print 'Target ' + det_id + ' was not found in the VLSS-catalogue'
	    #except:
		#print 'Target ' + det_id + ' was not found in the VLSS-catalogue, there was also no junk file'  
	    #else:
		#print 'Target ' + det_id + ' was not found in the VLSS-catalogue, there was also no junk file'  
	      
	    #if os.path.getsize(os.path.join(os.path.abspath("."), 'ClusterRegionImages/TGSS/TGSS-' + det_id + '.fits') ) < 10000:
		#os.remove(os.path.join(os.path.abspath("."), 'ClusterRegionImages/TGSS/TGSS-' + det_id + '.fits'))
		#print 'Target ' + det_id + ' was not found in the TGSS-catalogue'
		
	    #print >> file, "%6.0f"%det_id, "%12s"%ra_hms, "%13s"%dec_dms, "%6.2f"%z, "%8.2f"%L, "%6.2f"%scale, "%8.2f"%onempc
    return
    
Get_fits(infile)