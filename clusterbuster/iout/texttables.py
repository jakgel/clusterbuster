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
    """ Output to start a MCMC based analysis      
        Work in progress
    
    """
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
            unline += entry[0](protoobj).unit_string()
             
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
            """ Formats the line: For each object the entry ( a lambda variable) of the list is looked for 
            line += entry[1] %  entry[0](obj)  is the simple version, but because I allow for formating of several value sin one column, I
            
            """
            if isinstance(entry[0], list):
                """ Multiple entries """
                line += entry[1] % tuple([e(obj) for e in  entry[0]])  #else: mhhh
            elif len(entry)>3:
                """ One single entry """
                line += entry[1] %  entry[0](obj)  #else:
            else:
                """ Measurand """
                line += entry[1] %  entry[0](obj)()  #else:
            if ii < len(dictionary) - 1: line += delimiter       
        line += '%s \n' % (ender)
        columns.append(line)
        
    return "".join(columns)



def create_table(objectlist, dictionary, caption='nocap', outer=False, longtab=False):
    """ A function to create a late table based on an object list and an dictionary of values """
    
    header, ender = create_table_frame(objectlist[0], caption, dictionary, outer=outer, longtab=outer)
    columns       = create_table_columns(objectlist, dictionary)

    print(type(header), type(columns), type(ender))
    print(columns)
    return header + columns + ender


def RList2table_paper(location, survey, longtab=False):
    
    survey.FilterCluster()
    RList     = survey.fetch_totalRelics()
    RList.sort(key= iom.Object_natural_keys ) 
    
    
#'', r/ö/, .label, .unit, rrrr}  
#lambda x: cbclass.measurand( x.R200/1000       , '$R_{200}$', un = 'Mpc' )  
    dictionary = [ [lambda x: x.name.replace('_', ' '), '%25s', 'l', 'Identifier', ''],
                   [lambda x: x.RA , '%5.2f', 'r'],
                   [lambda x: x.Dec, '%+7.2f', 'r'],
                   #[lambda x: x.Mach, '%.1f', 'r'],
                   [lambda x: x.alpha, '%.2f', 'r'],
                   [[lambda x: x.flux(), lambda x: x.flux.std[0]], '$%7.1f\pm%5.1f$', 'r', '$S_{1.4}$', '[mJy]'],
                   [[lambda x: np.log10(x.P_rest()), lambda x: np.log10((x.P_rest()+x.P_rest.std[0])/x.P_rest()), lambda x:np.log10((x.P_rest()-x.P_rest.std[0])/x.P_rest())], '$%5.2f^{+%4.2f}_{%4.2f}$', 'r', 'log$_{10}(P_{1.4})$', '$\mathrm{[W/Hz^{-1}]}$'],
                   [lambda x: x.LAS, '%5.2f', 'r'], #add error if you like
                   [lambda x: x.LLS, '%5.0f', 'r'], #add error if you like
                   #[Omega not needed],
                   [lambda x: x.iner_rat(), '%.2f', 'r', '$\lambda_2/\lambda_1$', ''],
                   [lambda x: x.Dproj_pix, '%.0f', 'r'],
                   [lambda x: x.theta_rel(), '%.1f', 'r', '$\phi$', '[deg]']
                 ]
    caption = 'Radio relic emission islands identified in the %s images'  % (survey.name)
    table = create_table(RList, dictionary, caption=caption)
            
    # WRITE COMMON FILE
    mf = open(location,"w")
    mf.write(table)
    mf.close()


def GClList2table_paper(location, survey, shortlist=None, longtab=False):   
    

    outer = False
    n_clusters=0
    n_compact_sources = 0
  
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
    """ Start with relic cluster """
    mf = open(location,"w")
    mf.write(head)


    for m in survey.GCls:
        
        
        # WRITE COMMON FILE
        if m.getstatusstring()[0] is not None:
          
            # PS$^1$  &
            # \'\'/Mpc  &
            #   %4.1f &      1000./(m.cosmoPS*60)
        	
    #        if m.P_rest() == 0:
    #              print '!!!!!!!! m.Prest == 0', m.name,  m.P_rest
            string = "%25s & $%.3f$ &  %.1f & %s & $%5.1f$ &"   % (m.name, m.z.value, m.M200.value/1e14, m.getstatusstring()[1], m.flux_lit.value)  +   '%s' % (m.gettypestring())  + '& %s - %s\\\\' % (findshort(m.Lx.ref.ident, shortlist) , findshort(m.flux_lit.ref.ident, shortlist) )
            mf.write(string + '\n')
            n_clusters += 1

    print('n_clusters:',  n_clusters)

    mf.write( foot )
    mf.close()
    

    """ Write relic clusters coordinates """
    mf = open(location + "centre", "w")

    head_mod = """\\begin{table*}
                \\begin{center}
                \\caption{Radio relic hosting clusters centre coordinates}
                \\scriptsize""" * int(outer) + """ 
                \\begin{tabular}{ l  r  r  c  c  c  c}
                \hline\hline
                \label{tab:NVSSrelics}
                Cluster               &  RA  &  Dec    &  Method   &  \multicolumn{2}{c}{$\\Delta_\\mathrm{simbad}$} & References\\\\ 
                                      &  [deg]  &  [deg]  &           &       [\SI{}{\\arcsecond}]       &    $[R_{200}]$    &    \\\\ 
                    (1)               & (2)  &  (3)  &        (4)  &          (5)  &  (6)    &  (7)  \\\\\hline\hline     """

    foot_mod = '\hline\hline\n\\end{tabular}\n' + int(outer) * '\\end{center}\n\\label{tab:NVSSrelics}\n\\end{table*}'
    mf.write(head_mod)
    for m in survey.GCls:

        # WRITE COMMON FILE
        if m.getstatusstring()[0] is not None:
            # PS$^1$  &
            # \'\'/Mpc  &
            #   %4.1f &      1000./(m.cosmoPS*60)

            #        if m.P_rest() == 0:
            #              print '!!!!!!!! m.Prest == 0', m.name,  m.P_rest
            string = "%25s &  %.3f & %.3f &" % (m.name, m.RA.value, m.Dec.value) + 'Planck & 0.01  & 0.1  & cite \\\\'
            mf.write(string + '\n')

    mf.write(foot_mod)
    mf.close()


    """ Write subtracted sources """
    mf = open(location + "compacts", "w")

    head_mod = """\\begin{table*}
                \\begin{center}
                \\scriptsize""" * int(outer) + """ 
                \\begin{tabular}{ l  r  r  r  c  c  c  c}
                \hline\hline
                \label{tab:NVSSrelics}
                Cluster               &  RA     &  Dec    &  flux  &  type &  $\\Theta_\\mathrm{major}$ & $\\Theta_\\mathrm{minor}$ & $\\theta$ \\\\ 
                                      &  [deg]  &  [deg]  &  [mJ]  &              &      [ \SI{}{\\arcminute}]  &            [\SI{}{\\arcminute}]             &    [deg]    \\\\ 
                    (1)               & (2)     &  (3)    &   (4)  &     (5)      &       (6)                  &  (7)                      &    (8)     \\\\\hline\hline
                """

    foot_mod = '\hline\hline\n\\end{tabular}\n' + int(outer) * '\\end{center}\n\\label{tab:NVSS_compactsources}\n\\end{table*}'
    mf.write(head_mod)
    for m in survey.GCls:

        # WRITE COMMON FILE
        if m.getstatusstring()[0] is not None:

            for compact in m.compacts:
                values = iom.J2000ToCoordinate(compact['dir'])
                values = [float(value) for value in values]
                if values[0] > 0:
                    RA = (values[0]+values[1]/60+values[2]/3600)*15
                else:
                    RA = (values[0]-values[1]/60-values[2]/3600)*15
                Dec = values[3]+values[4]/60+values[5]/3600
                string = "%25s &  %.3f & %.3f &" % (m.name, RA, Dec)
                if compact['shape'] == "Point":
                    string += "%.1f & %s  &  &   &  \\\\" % (float(compact['flux']), compact['shape'])
                else:
                    string += "%.1f & %s & %.2f & %.2f  & %.0f \\\\" % (float(compact['flux']), "Extended",
                                                                            float(compact['majoraxis']), float(compact['minoraxis']),
                                                                            float(compact['theta']))

                mf.write(string + '\n')
                n_compact_sources += 1

    print('n_compact_sources:', n_compact_sources)

    mf.write(foot_mod)
    mf.close()


    """ Do phoenix clusters """
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
            string = "%25s & $%.3f$ &  %.1f & %s & $%5.1f$ &"   % (m.name, m.z.value, m.M200.value/1e14, m.getstatusstring()[1], m.flux_lit.value)  +   '%s' % (m.gettypestring())  + '& %s - %s\\\\' % (findshort(m.Lx.ref.ident, shortlist) , findshort(m.flux_lit.ref.ident, shortlist) )
            mf.write(string + '\n')
    
    mf.write( foot )
    mf.close()




   
def LOFARwiki(survey):
    """ A functionality created to embed images into an internal wiki """
    
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