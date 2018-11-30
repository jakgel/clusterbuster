#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:54:38 2017

@author: jakobg
"""

#!/usr/bin/env python

# Some classes for the relic lsit
# ...
from __future__ import division,print_function

import csv
import numpy as np
import clusterbuster.surveyclasses as cbclass
import clusterbuster.dbclasses     as dbclass
import clusterbuster.iout.misc     as iom
import os



def findshort(phrase, shortlist):
    
    if shortlist is not None:
        
        if phrase is '':
            return '\phantom{AAA00}'
        
        with open(shortlist, 'rb') as csvfile:        
            shortlist = csv.DictReader(csvfile, delimiter=',', quotechar='"')

            for CL in shortlist:
                if phrase==CL['Long']: return CL['Short']
 
    return phrase

def ClusterVol_lum(surveys, location):  
    ''' Output to start a MCMC based analysis      
        Work in progress
    
    '''
    GCllist = []
    header  = 'SurveyID, PowerCluster, M200, Simu?'
    for ii,survey in enumerate(surveys):
        print('ClusterVol_lum:', survey.name)
        for GCl in survey.GCls:
            
            if survey.Rmodel.simu:
                #V_vol:   
                if 'MUSIC2' in survey.name:
                    array = [ii, GCl.Prest_vol, GCl.M200.value, survey.Rmodel.simu, len(survey.GCls)] #[ii, GCl.Prest_vol.val, survey.Rmodel.simu]
                    GCllist.append(array)
                else:
                    array = [ii, GCl.Prest_vol.value, GCl.M200.value, survey.Rmodel.simu, len(survey.GCls)] #[ii, GCl.Prest_vol.val, survey.Rmodel.simu]
                    GCllist.append(array)       
            else:
                array = [ii, GCl.P_rest.value, GCl.M200.value,  survey.Rmodel.simu, len(survey.GCls)]
                GCllist.append(array)

    print("Save results in %s%s.csv" % (surveys[0].outfolder, location))
    np.savetxt("%s%s.csv" % (surveys[0].outfolder, location), GCllist, delimiter=",",  fmt='%.3e', header=header) 
            
    return 0
        
        




def create_table_frame(protoobj, caption, dictionary, delimiter='&', ender='\\', longtab=False, outer=True):
    
    # l  r  r  r  r  r  r  r  r  r  r}
    
    
    
    
    tabline  = ''
    idline   = ''
    unline   = ''
    sizeline = '\\scriptsize\n'
    labeline = "" #"'\label{tab:%s}\n' % (label)
    capline  = '\caption{%s} \n' % (caption)
    delime   = '\hline\hline\n'
    
    for ii,entry in enumerate(dictionary):
        tabline += ' %s ' % entry[2]
        if  not isinstance(entry[0], list) and isinstance(entry[0](protoobj), dbclass.measurand):
            idline += entry[0](protoobj).label
            unline += entry[0](protoobj).unit
             
        else:
             idline += entry[3]
             unline += entry[4]
        if ii < len(dictionary) - 1: 
            idline += delimiter 
            unline += delimiter 
     
    idline += ender + ender + '\n'
    unline += ender + ender #+ '\hline\hline'
    tabline += '} \n'
    
    
    
    
    
    if longtab:
        tabline  = '\\begin{longtable}{'+tabline
        capline += """\\\\"""
        head = """\\begin{landscape}\n\\begin{center}""" + sizeline + tabline + capline + delime + labeline +idline +unline + """
                            \\endfirsthead
                            \multicolumn{11}{l}{\\tablename\ \\thetable\ -- \\textit{Continued from previous page}} \\\\
                            \hline""" + idline +unline + """
                            \\endhead
                            \hline\hline \multicolumn{11}{r}{\\textit{Continued on next page}} \\\\
                            \\endfoot""" + delime + '\\endlastfoot'
        foot = '\hline\n\\end{longtable}\n\\end{center}\n\\end{landscape}'  
    else:
        
        tabline =  '\\begin{tabular}{' + tabline
        
        
        head =  ("""\\begin{table*}\n\\begin{center}""" + capline + sizeline) * int(outer) + tabline + delime +labeline + idline + unline + delime
        foot = '\hline\hline\n\\end{tabular}\n'+int(outer)*'\\end{center}\n\\end{table*}'
                                   

        
                   
#    lines = [taline,delime, idline,unlune]
#    head  = "".join(lines)    
    
    return head, foot




def create_table_columns(objectlist, dictionary, delimiter='&', ender='\\\\'):
    
    columns = []
    for obj in objectlist:
        line= ''
        for ii,entry in enumerate(dictionary):
            ''' Formats the line: For each object the entry ( a lambda variable) of the list is looked for 
            line += entry[1] %  entry[0](obj)  is the simple version, but because I allow for formating of several value sin one column, I
            
            '''
            if isinstance(entry[0], list):
                ''' Multiple entries '''
                line += entry[1] % tuple([e(obj) for e in  entry[0]])  #else: mhhh
            elif len(entry)>3:
                ''' One single entry '''
                line += entry[1] %  entry[0](obj)  #else:
            else:
                ''' Measurand '''
                line += entry[1] %  entry[0](obj)()  #else:
            if ii < len(dictionary) - 1: line += delimiter       
        line += '%s \n' % (ender)
        columns.append(line)
        
    return "".join(columns)



def create_table(objectlist, dictionary, caption='nocap', outer=False, longtab=False):
    ''' A function to create a late table based on an object list and an dictionary of values '''
    
    header, ender = create_table_frame(objectlist[0], caption, dictionary, outer=outer, longtab=outer)
    columns       = create_table_columns(objectlist, dictionary)

    print(type(header), type(columns), type(ender))
    print(columns)
    return header + columns + ender


def RList2table_paper_new(location, survey, longtab=False):   
    
    survey.FilterCluster(zborder=0.05,ztype='>')
    RList     = survey.fetch_totalRelics()
    RList.sort(key= iom.Object_natural_keys ) 
    
    
#'', r/ö/, .label, .unit, rrrr}  
#lambda x: cbclass.measurand( x.R200/1000       , '$R_{200}$', un = 'Mpc' )  
    dictionary = [ [lambda x: x.name.replace('_',' '), '%25s', 'l' , 'Identifier', ''],
                   [lambda x: x.RA , '%5.2f' , 'r'],
                   [lambda x: x.Dec, '%+7.2f', 'r'],
                   [lambda x: x.Mach, '%.1f', 'r'],
                   [[lambda x: x.flux(), lambda x: x.flux.std[0]], '$%7.1f\pm%5.1f$' , 'r', '$S_{1.4}$', '[mJy]'],
                   [[lambda x: np.log10(x.P_rest()), lambda x: np.log10((x.P_rest()+x.P_rest.std[0])/x.P_rest()), lambda x:np.log10((x.P_rest()-x.P_rest.std[0])/x.P_rest())], '$%5.2f^{+%4.2f}_{%4.2f}$', 'r', 'log$_{10}(P_{1.4})$' , '[W/Hz$^{-1}$]'],
                   [lambda x: x.LAS, '%5.2f', 'r'], #add error if you like
                   [lambda x: x.LLS, '%5.0f', 'r'], #add error if you like
                   #[Omega not needed],
                   [lambda x: x.iner_rat(), '%.2f', 'r', '$\lambda_2/\lambda_1$', '[]'],
                   [lambda x: x.Dproj_pix, '%.0f', 'r'],
                   [lambda x: x.theta_rel(), '%.1f', 'r', '$\phi$', '[deg]']
                 ]
    caption = 'Radio relic emission islands identified in the %s images'  % (survey.name)
    table = create_table(RList, dictionary, caption=caption)
            
    # WRITE COMMON FILE
    mf = open(location,"w")
    mf.write(table)
    mf.close()
  
    
    
    
def RList2table_paper(location, survey, longtab=False):   

    RList    = survey.fetch_totalRelics(zborder=0.05)
    # To sort the list in place..., sort for Cluster ID
    RList.sort(key= iom.Object_natural_keys ) 
    outer = False 
  
    if longtab:
        head = """\\begin{landscape}
    \\begin{center}
    \\scriptsize
    \\begin{longtable}{ l  r  r  r  r  r  r  r  r  r  r}
    % \centering
    \caption{Radio relic emission islands identified in the NVSS images}\\\\ \hline\hline
    \label{tab:NVSSrelics}
    Identifier                &  RA    &  Dec   &  $S_{1.4}$  &   log$_{10}(P_{1.4})$        &   LAS    &  LLS   & $\Omega$       &  $\\lambda_2/\\lambda_1$   & $D_\mathrm{proj,pix}$  &  $\phi$ \\\\ 
                              &  [deg] & [deg]  & [mJy]       &   [W/Hz$^{-1}$]   & [arcmin] &  [kpc]   &  [$\Omega_\mathrm{beam}$]  &   []    &  [kpc]     &  [deg]     \\\\ \hline
    \\endfirsthead
    \multicolumn{11}{l}{\\tablename\ \\thetable\ -- \\textit{Continued from previous page}} \\\\
    \hline
    Identifier                &  RA    &  Dec   &    $S_{1.4}$    &   log$_{10}(P_{1.4})$     &   LAS & LLS   &  $\Omega$      &  $\\lambda_2/\\lambda_1$   & $D_\mathrm{proj,pix}$  &  $\phi$ \\\\ 
                              &  [deg] & [deg]  & [mJy]       &   [W/Hz$^{-1}$]   & [arcmin] &  [kpc]   &  [$\Omega_\mathrm{beam}$]  &   []    &  [kpc]     &  [deg]     \\\\ \hline
    \\endhead
    \hline\hline \multicolumn{11}{r}{\\textit{Continued on next page}} \\\\
    \\endfoot
    \hline\hline
    \\endlastfoot"""
        foot = '\hline\n\\end{longtable}\n\\end{center}\n\\end{landscape}'  
    else:
        head =    """\\begin{table*}
    \\begin{center}
    \caption{Radio relic emission islands identified in the NVSS images}
    \\scriptsize""" * int(outer) + """ 
    \\begin{tabular}{ l  r  r  r  r  r  r  r  r  r  r}
    \hline\hline
    \label{tab:NVSSrelics}
    Identifier                &  RA    &  Dec   &    $S_{1.4}$ &   log$_{10}(P_{1.4})$     &   LAS  &  LLS   & $\Omega$       &  $\\lambda_2/\\lambda_1$   & $D_\mathrm{proj,pix}$  &  $\phi$ \\\\ 
                              &  [deg] & [deg]  & [mJy]       &   [W/Hz$^{-1}$]   & [arcmin] &  [kpc]   &  [$\Omega_\mathrm{beam}$]  &   []    &  [kpc]     &  [deg]     \\\\ \hline\hline 
    """
        foot = '\hline\hline\n\\end{tabular}\n'+int(outer)*'\\end{center}\n\\label{tab:NVSSrelics}\n\\end{table*}'  
    
    
    # WRITE COMMON FILE
    mf = open(location,"w")
    mf.write(head)
  
  
#  formatingDic = {  name, Identifier, 'deg', lambda x: x.namestring  "%25s & %5.2f & %+7.2f  & $%7.1f\pm%5.1f$ &",
#                    name,     .label,  .un, lambda x: x.namestring  "%25s & %5.2f & %+7.2f  & $%7.1f\pm%5.1f$ &", }

# PS$^1$  &
# \'\'/Mpc  &
#   %4.1f &      1000./(m.cosmoPS*60)
    name      = ''
    clname    = ''
    for ii,m in enumerate(RList):
        if m.region.rtype.classi > 0:
      
            if name ==  m.name.replace('_',' '):
                namestring = " " * len(name)
            else:
                namestring = m.name.replace('_',' ')
                name       = namestring
#    if m.GCl.name ==  clname:
##	   lenname = 'A' * 99
##	   for alias in aliases:
##	      clname_a     =  m.GCl.name.replace('_',' ')
##	      clname_b     =  clname_a.replace(alias[0], alias[1]) 
##	      if len(clname_b) < len(lenname):
##		 lenname  = clname_b
#      namestring = namestring.replace( clname.replace('_',' ') , '\\phantom{%s}' % ( lenname )  ) # phantom{} if name repeats ...
#	  
#           
#      for alias in aliases:
#          namestring = namestring.replace(alias[0], alias[1]) 
#	
#      clname = m.GCl.name 
	
	
        if m.P_rest() == 0:
            print('!!!!!!!! m.Prest == 0', m.name,  m.P_rest)
            string = "%25s & %5.2f & %+7.2f  & $%7.1f\pm%5.1f$ &"   % (namestring, m.RA(), m.Dec(), m.flux(), m.flux.std[0])  \
                     +  "$%5.2f^{+%4.2f}_{%4.2f}$ &  $%5.2f$&" % (np.log10(m.P_rest()), np.log10((m.P_rest()+m.P_rest.std[0]/m.P_rest())), np.log10((m.P_rest()-m.P_rest.std[1])/m.P_rest()), m.LAS()) \
                     +  "$%5.0f\pm%5.0f$ & $%5.1f\pm%3.1f$ &" %( m.LLS(), m.LLS.std[0], m.area()/m.dinfo.Abeam[1], m.area.std[0]/m.dinfo.Abeam[1] )   \
                     +  "%5.2f  &  %4.0f  &  %5.1f \\\\" %(m.iner_rat(), m.Dproj_pix(), m.theta_rel())
      
        #string = "%25s  %5.1f \\" %(namestring, m.theta_rel[1]) 
      
        mf.write(string + '\n')
    
    mf.write( foot )
    #mf.write('\hline\n\\end{tabular}\n\\end{table}')
    mf.close()


def GClList2table_paper(location, survey, shortlist=None, longtab=False):   
    

    outer = False 
  
    if longtab:
        head = """\\begin{landscape}
\\begin{center}
\\scriptsize
\\begin{longtable}{ l  r  r  r  r  c  c}
% \centering
\caption{Radio relic hosting clusters}\\\\ \hline\hline
\label{tab:NVSSrelics}
Cluster               &  z   &  $M_{200}$          &  $F_\mathrm{NVSS}$   &   $F_\mathrm{1.4, lit}$      &  Diff. Emi.   &  References\\\\ 
                      & []   &  [$10^{14} M_\odot$]  &       [mJy]          &          [mJy]               &  RHR*         &     mass $\\mid$ relics  \\\\ 
    (1)               & (2)  &  (3)                &        (4)           &          (5)                 &  (6)          &     (7)    \\\\\hline\hline         
\\endfirsthead
\multicolumn{11}{l}{\\tablename\ \\thetable\ -- \\textit{Continued from previous page}} \\\\
\hline
Cluster               &  z   &  $M_{200}$          &  $F_\mathrm{NVSS}$   &   $F_\mathrm{1.4, lit}$      &  Diff. Emi.   &  References\\\\ 
                      & []   &  [$10^{14} M_\odot$]  &       [mJy]            &          [mJy]               &  RHR*         &     mass $\\mid$ relics  \\\\ 
    (1)               & (2)  &  (3)                &        (4)           &          (5)                 &  (6)          &     (7)    \\\\\hline\hline         
\\endhead
\hline\hline \multicolumn{11}{r}{\\textit{Continued on next page}} \\\\
\\endfoot
\hline\hline
\\endlastfoot"""
        foot = '\hline\n\\end{longtable}\n\\end{center}\n\\end{landscape}'  
    else:
        head =    """\\begin{table*}
\\begin{center}
\caption{Radio relic hosting clusters}
\\scriptsize""" * int(outer) + """ 
\\begin{tabular}{ l  c  r  r  r  c  c}
\hline\hline
\label{tab:NVSSrelics}
Cluster               &  z   &  $M_{200}$          &  $F_\mathrm{NVSS}$   &   $F_\mathrm{1.4, lit}$      &  Diff. Emi.   &  References\\\\ 
                      & []   &  [$10^{14} M_\odot$]  &       [mJy]            &          [mJy]               &  RHR*         &     mass $\\mid$ relics  \\\\ 
    (1)               & (2)  &  (3)                &        (4)           &          (5)                 &  (6)          &     (7)    \\\\\hline\hline                  
                      
"""
        foot = '\hline\hline\n\\end{tabular}\n'+int(outer)*'\\end{center}\n\\label{tab:NVSSrelics}\n\\end{table*}'  
    


    
    
    #RList.sort(key= iom.Object_natural_keys ) 
    ''' Start with relic cluster '''
    mf = open(location,"w")
    mf.write(head)

    n_clusters=0
    for m in survey.GCls:
        
        
        # WRITE COMMON FILE
        if m.getstatusstring()[0] is not None:
          
            # PS$^1$  &
            # \'\'/Mpc  &
            #   %4.1f &      1000./(m.cosmoPS*60)
        	
    #        if m.P_rest() == 0:
    #              print '!!!!!!!! m.Prest == 0', m.name,  m.P_rest
            string = "%25s & $%.3f$ &  %.2f & %s & $%5.1f$ &"   % (m.name, m.z.value, m.M200.value/1e14, m.getstatusstring()[1], m.flux_lit.value)  +   '%s' % (m.gettypestring())  + '& %s - %s\\\\' % (findshort(m.Lx.ref.ident, shortlist) , findshort(m.flux_lit.ref.ident, shortlist) )
            mf.write(string + '\n')
            n_clusters += 1

    print('n_clusters:',  n_clusters )

    mf.write( foot )
    mf.close()
    
    
    ''' Do phoenix clusters '''
    mf = open(location+'phoenixes',"w")
    mf.write(head)
    for m in survey.GCls:
        
#        if m.flux() == 0 or m.z < 0.05:
        if m.getstatusstring()[0] is None:
                
            # PS$^1$  &
            # \'\'/Mpc  &
            #   %4.1f &      1000./(m.cosmoPS*60)
        	
    #        if m.P_rest() == 0:
    #              print '!!!!!!!! m.Prest == 0', m.name,  m.P_rest
            string = "%25s & $%.3f$ &  %.2f & %s & $%5.1f$ &"   % (m.name, m.z.value, m.M200.value/1e14, m.getstatusstring()[1], m.flux_lit.value)  +   '%s' % (m.gettypestring())  + '& %s - %s\\\\' % (findshort(m.Lx.ref.ident, shortlist) , findshort(m.flux_lit.ref.ident, shortlist) )
            mf.write(string + '\n')
    
    mf.write( foot )
    mf.close()



   
def LOFARwiki(survey):
    ''' A functionality created to embed images into an internal wiki '''
    
    # Create LOFAR-wiki links
    string = ''
    for ii, GCl in enumerate(survey.GCls):
        if ii % 3  == 0:
            string += '\n'      
        string += '<<EmbedObject(NVSS-%s.pdf,width=500.0,height=500,menu = False)>>\n' % (GCl.name)
        
    return string
    
#=== Old===#
#==========#
def RList2nparray_csv(RList, AddInfo):  
    
    
  nparray =  np.empty( (19, len(RList)), dtype=str) 
    
  for ii,relic in enumerate(RList):
      nparray[:,ii] = np.asarray( [relic.GCl.name, relic.region.name, relic.GCl.z, relic.RA, relic.Dec, relic.region.rtype, relic.dinfo.limit, relic.flux, np.log10(relic.P), np.log10(relic.Prest), relic.Dproj[0], relic.Dproj[1], relic.LLS, relic.LAS, 0, 0, 0, 0, relic.area] )

  header     = "ClusterID; relic; z; RA; Dec; rtype; detlimit; flux; log10(Pobsframe[W/Hz/m^2]); log10(P1400MHz); Distance(arcsec); Distance(kpc); LLS; LAS; a (kpc); b (kpc); a(arcsec); b(arcsec); Area(kpc^2); Area(arcsec^2); angle(relic); relative angle; alpha; alpha?"
  
  return nparray, header


def RList2dat_Seb_relic(location, GCls):   
    
   
  # To sort the list in place..., sort for Cluster ID
  GCls.sort(key=lambda GCl: GCl.name, reverse=False) 
  mf        = open(location,"w")
  mf.write("""#              name                   z       RA       DEC       L_x    flux   flux_err     P_14     P_14_err   LLS     LLS_err  area_amin2     D_proj   Class flux_lit  LLS_lit  L_x_lit area_kpc2   alpha alpha_err alpha_FLAG   scale    a_GCL a_elo_cnt a_elo_iner a_elo_flux a_rel_cnt a_rel_iner a_rel_flux class_FLAG   ecc  area_amin2_err area_100kpc_err
    #                                                          10^44erg/s    mJy       mJy   10^24W/Hz  10^24W/Hz   Mpc      Mpc       arcmin2       kpc                            10^24W/Hz    100kpc2                               '/Mpc     deg      deg        deg        deg       deg        deg        deg                 arcmin2    100kp2
    #              S50                   f8       f8       f8        f8       f8       f8       f8         f8       f8       f8           f8          f8       i4      f8       f8       f8       f8        f8       f8          i4     f8        f8       f8         f8         f8        f8         f8         f8          i4     f8     f8         f8
    #               0                     1        2        3         4        5        6        7          8        9       10           11          12       13      14       15       16       17        18       19          20     21        22       23         24         25        26         27         28          29     30      31        32
    #\n""")
      #mf.write("#id [0]   &  z  [1]     &    RA  [2]  &  DEC [3]   & $L_x[0.1-2.3] (10^44 erg/s) [4] $  & Flux (mJy) [5] & Flux_err [6] & $P_{1.4, rest} (10^24 W/Hz) [7] $ & P_{1.4, rest}_error (10^24 W/Hz)  [8] &  LLS (Mpc) [9] & LSS_error [10] & Area (arcmin^2) [11]&   D_proj (kpc) [12]& Classification [13]& Flux_lit  (mJy) [14] & LLS_lit [15] & Lx_lit[0.1-2.3] [16] & area (10^4 kpc^2) [17]&  alpha  [18] & alpha_err  [19]  & alpha_FLAG  [20] &   '/Mpc [21] & angle_GCl [22] & angle_elo (cnt) [23] & angle_elo (iner) [24] & angle_elo (flux) [25] & angle_rel (cnt) [26] & angle_rel (iner) [27] & angle_rel (flux) [28] & Class_FLAG [29] & iner_rat (iner) [30]\n")  
      
  for GCl in GCls:
      relics = GCl.relics
      relics.sort(key=lambda relics: relics.region.name, reverse=False)
      for m in relics:             
        string =  "%27s %8.3f %8.3f %8.3f %8.2f" %(m.name, GCl.z, m.RA, m.Dec, GCl.Lx) +  " %8.2f %8.2f %8.4f %8.4f" %(m.flux, m.flux_err, m.Prest/1.e24, m.Prest_err/1.e24) +  " %8.4f %8.4f %8.3f %8.1f" % (1e-3*m.LLS, 1e-3*m.LLS_err, m.area, m.Dproj) + " %2i %8.2f %8.2f %8.2f" % (m.region.rtype.classi, m.flux_lit, m.LLS_lit, GCl.Lx_lit) + " %8.2f  %8.2f %8.2f %i" % (m.area100kpc, m.region.alpha, m.region.alpha_err, m.region.alphaFLAG) + " %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %3i %8.4f" %(1e3/(GCl.cosmoPS*60), m.theta_GCl,  m.theta_elong[0], m.theta_elong[1], m.theta_elong[2], m.theta_rel[0], m.theta_rel[1], m.theta_rel[2], GCl.ClassFlag, m.ecc) + " %7.4f %7.4f" %(m.area_err, m.area100kpc_err)
        mf.write(string + '\n') 
    
  mf.close()  
 
def RList2dat_MUSIC2_relic(location, RList, Bmodel, append=False):   
    
   
  # To sort the list in place..., sort for Cluster ID
  (relics,eff,mockobs) = RList
  
  header = """#              name                   z       RA       DEC       L_x    flux   flux_err     P_14     P_14_err   LLS     LLS_err  area_amin2     D_proj   Class flux_lit  LLS_lit  L_x_lit area_kpc2   alpha alpha_err alpha_FLAG   scale    a_GCL a_elo_cnt a_elo_iner a_elo_flux a_rel_cnt a_rel_iner a_rel_flux class_FLAG   ecc  area_amin2_err area_100kpc_err   eff  Bmodel.idgrid[0]  Bmodel.idgrid[1]  mockobs.theta  mockobs.psi  mockobs.phi  mockobs.snap  mockobs.clid 
#                                                          10^44erg/s    mJy       mJy   10^24W/Hz  10^24W/Hz   Mpc      Mpc       arcmin2       kpc                            10^24W/Hz    100kpc2                               '/Mpc     deg      deg        deg        deg       deg        deg        deg                 arcmin2    100kp2
#              S50                   f8       f8       f8        f8       f9       f9       f8         f8       f8       f8           f8          f8       i4      f8       f8       f8      f10        f8       f8          i4     f8        f8       f8         f8         f8        f8         f8         f8          i4     f8     f9         f9     e8     i4     i4    f8    f7   f7    i3    i5
#               0                     1        2        3         4        5        6        7          8        9       10           11          12       13      14       15       16       17        18       19          20     21        22       23         24         25        26         27         28          29     30      31        32     33     34     35    36    37   38    39    40
#\n"""
  
  if append:
     try:
        os.stat( location )
        mf        = open(location,"a")
     except:
        mf        = open(location,"w")
        mf.write(header)
  else:
        mf        = open(location,"w")
        mf.write(header)
        
  for ii,m in enumerate(relics):
    #print  (m.GCl.cosmoPS*60/100)**2               
    string =  "%27s %8.3f %8.3f %8.3f %8.2f" %(m.name, m.GCl.z, m.RA, m.Dec, m.GCl.Lx) +  " %9.2f %9.2f %8.4f %8.4f" %(m.flux, m.flux_err, m.Prest/1.e24, m.Prest_err/1.e24) +  " %8.4f %8.4f %8.3f %8.1f" % (1e-3*m.LLS, 1e-3*m.LLS_err, m.area, m.Dproj[1]) + " %2i %8.2f %8.2f %8.2f" % (m.region.rtype.classi, m.flux_lit, m.LLS_lit, m.GCl.Lx_lit) + " %10.2f  %8.2f %8.2f %i" % (m.area100kpc, m.region.alpha, m.region.alpha_err, m.region.alphaFLAG) + " %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %3i %8.4f" %(1e3/(m.GCl.cosmoPS*60), m.theta_GCl,  m.theta_elong[0], m.theta_elong[1], m.theta_elong[2], m.theta_rel[0], m.theta_rel[1], m.theta_rel[2], m.GCl.ClassFlag, m.ecc) + " %9.4f %9.4f" %(m.area_err, m.area100kpc_err) + "%11.4e %4i %4i" %(eff, Bmodel.idgrid[0], Bmodel.idgrid[1]) + "%8.4f %7.4f %7.4f %3i %5i" %(mockobs.theta, mockobs.psi, mockobs.phi, mockobs.snap, mockobs.clid) 
    mf.write(string + '\n') 

  mf.close()  
  
 
def RList2dat_Seb_cluster(location, ClList):   
    
  mf        = open(location,"w")
  mf.write("""#                   name          z       RA      DEC       L_x      flux    flux_err     P_14    P_14_err   area_amin2  area_kpc2   alpha  alpha_err alpha_FLAG  D_proj   flux_lit   Class  firmClassi?  area_amin2_err area_100kpc_err
#                                                       10^44 erg/s   mJy       mJy    10^24 W/Hz 10^24 W/Hz   arcmin2     kpc2                                    kpc         mJy                        arcmin2 100kpc2
#                    S50         f8      f8       f8        f8        f8        f8        f8         f8          f8         f8         f8        f8       S5        f8         f8      5s       5s       f8      f8
#                     0           1       2        3         4         5         6         7          8           9         10         11        12       13        14         15      16       17       18      19
#\n""")
  #mf.write("#id [0]   &  z  [1]     &    RA  [2]  &  DEC [3]   & $L_x[0.1-2.3] (10^44 erg/s) [4] $  & Flux (mJy) [5] & Flux_err [6] & $P_{1.4, rest} (10^24 W/Hz) [7] $ & P_{1.4, rest}_error (10^24 W/Hz)  [8] &   Area (arcmin^2) [9] &   Area (100kpc^2) [10] &  alpha  [11] & alpha_err  [12]  & alpha_FLAG  [13] &  D_proj (kpc) [14]& F_lit (mJy) [15] & Class [16] & firm classi? [17]\n")  
  for ii,cl in enumerate(ClList):
    string = "%27s %8.3f %8.3f %8.3f"  % (cl.name, cl.z, cl.RA, cl.Dec) + " %8.3f %8.3f %8.3f %9.5f %8.5f"  %(cl.Lx, cl.flux, cl.flux_err, cl.Prest/1.e24, cl.Prest_err/1.e24)  + " %8.3f %8.3f %8.2f %8.2f %6s" % (cl.area, cl.area100kpc, cl.region.alpha, cl.region.alpha_err, cl.region.alphaFLAG) + " %8.3f %6.1f %5s %5s" %(cl.Dproj[1], cl.flux_lit, cl.region.rtype.classi, cl.ClassFlag) + " %7.4f %7.4f" %(cl.area_err, cl.area100kpc_err)
    mf.write(string + '\n')
  mf.close() 
