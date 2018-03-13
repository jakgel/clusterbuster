import numpy                as np
import scipy.interpolate    as interpolate
import matplotlib.pyplot    as plt
import pyutil_my.MATHutil   as math
'''
    Start with e.g. InterpolateRadio2D(psiFile = '../Analysis_MUSIC2/Hoeft_radio/mach_psi_tablefine(10,3).txt',  inter=(10,6)) 
'''



# from  http://stackoverflow.com/questions/5328128/scipy-interpolation-of-large-matrix
def my_interp(X, Y, Z, x, y, spn=3):
    xs,ys = map(np.array,(x,y))
    z = np.zeros(xs.shape)
    for i,(x,y) in enumerate(zip(xs,ys)):
        # get the indices of the nearest x,y
        xi = np.argmin(np.abs(X[0,:]-x))
        yi = np.argmin(np.abs(Y[:,0]-y))
        xlo = max(xi-spn, 0)
        ylo = max(yi-spn, 0)
        xhi = min(xi+spn, X[0,:].size)
        yhi = min(yi+spn, Y[:,0].size)
        # make slices of X,Y,Z that are only a few items wide
        nX = X[xlo:xhi, ylo:yhi]
        nY = Y[xlo:xhi, ylo:yhi]
        nZ = Z[xlo:xhi, ylo:yhi]
        intp = interpolate.interp2d(nX, nY, nZ)
        z[i] = intp(x,y)[0]
    return z


# from here on: done by myself

def LoadFile_psi(psiFile):
    ''' Just gives the Mach number and Temeprature values '''
    #=== FILE A ===#
    # read first line .... split it and convert sstring to float science float('1.31E+01') or for a list:map(float, ['3.76E+00', '1.31E+01', '1.14E+01'])


    with open(psiFile, 'r') as f:
        first_line = f.readline()
        psi_x      = first_line.split()[2:]                # Splits into list without first two elements
    psi_x      = np.asarray( [float(i) for i in psi_x ]  ) # Converts strings to floats   # Converts strings to floats
    psi_y       = np.loadtxt(psiFile,skiprows=0)[:,0]
    
    return psi_x, psi_y


def InterpolateRadio2D(psiFile='../Analysis_MUSIC2/Hoeft_radio/mach_psi_table.txt', machFile='../Analysis_MUSIC2/Hoeft_radio/q_mach_machr_table.txt', saveplot='../Analysis_MUSIC2/Hoeft_radio/interpolated', psiFileNew = False, machFileNew = False,  inter=(10,3)):
# Currently the mach number is interpolated in an logarithmic space which is much sparser at lower mach numbers then anticipated  
# I suspect an double-exponential function for mach (both efficiency dependency stepsize)  
  
    # Note that the original grid given in 'Hoeft_radio/mach_psi_table.txt' is (quite) regular in log-loglog space, which makes it very simple to invoke an interpolation function!
    # Irregular data points would make it nececcairy to use functions like scipy.interpolate.griddata(points, values, (grid_x, grid_y), method='cubic')
  
    plot_old = False
    plot_new = False
    plot_PhD = True
    
    ##==== psiFile   for psi factor; machfile for mach-numbers conversion factors
    H_mach      = np.loadtxt(machFile,skiprows=0) 
    H_psi       = np.loadtxt(psiFile,skiprows=0)[:,1::] # you wont get the temperature values ... read them separetely
    psi_x,psi_y = LoadFile_psi(psiFile)
    
    psi_x      = np.log10( psi_x )             # converts to and log10        space
    psi_y       = np.log10(np.log10( psi_y ))  # converts to and log10(log10) space
    X, Y = np.meshgrid(psi_x, psi_y)
    Z = np.log10(H_psi)
    
    #interp_spline   =  interpolate.interp2d(x, y, Z) #, kind='cubic'
    interp_spline    =  interpolate.RectBivariateSpline(psi_y, psi_x, Z) #, bbox=[None, None, None, None], kx=3, ky=3, s=0
        
    xnew = np.arange(psi_x[0], psi_x[-1], (psi_x[-1]-psi_x[0])/(len(psi_x)*inter[0]) ) #np.arange(-4,  2, 4e-2) #
    ynew = np.arange(psi_y[0], psi_y[-1], (psi_y[-1]-psi_y[0])/(len(psi_y)*inter[1]) ) #np.arange(0.2, 3, 2e-2) #
    Znew = interp_spline(ynew, xnew )
    
    keV2K = 1.16e7 # Translates keV to Kelvin
    
    if plot_old:
        print '!!!'
        plt.plot( np.arange(0, len(psi_x), 1 ), psi_x )
        plt.plot( np.arange(0, len(psi_y), 1 ), psi_y )
        plt.savefig(saveplot + '_linearity.png')
        
        fig = plt.figure()
        
        ax1 = plt.subplot(121)
        ax1.pcolor( np.log10(keV2K) + psi_x, psi_y, Z)
        ax1.set_title("Sparsely sampled function")
        ax1.set_xlim([3.1, 9])
        ax1.set_ylim([psi_y[0], 0.5])
        ax1.set_xlabel('$\\mathrm{log_{10}(T)\\,[K]}$ ')
        ax1.set_ylabel('$\\mathrm{log_{10}(log_{10}(M))\\,[]}$')
    
    
        ax2 = plt.subplot(122)
        im2 = ax2.pcolor( np.log10(keV2K) + xnew, ynew, Znew)
        ax2.set_title("Interpolated function")
        ax2.set_xlim([3.1, 9])
        ax2.set_ylim([psi_y[0], 0.5])
        ax2.set_xlabel('$\\mathrm{log_{10}(T)\\,[K]}$ ')
        ax2.set_yticklabels([])
    
        mach  = [1.5,2.2,3.0,10.0]
        c     = [plt.cm.rainbow( (np.log10(np.log10(m))-ax1.get_ylim()[0])/abs(ax1.get_ylim()[1]-ax1.get_ylim()[0]) ) for m in mach]
        for ii,m in enumerate(mach):
            ax1.plot( [ax1.get_xlim()[0], ax1.get_xlim()[1]] , [np.log10(np.log10(m))]*2, '-', c=c[ii], lw=1.5, alpha=0.9 ) 
            ax2.plot( [ax2.get_xlim()[0], ax2.get_xlim()[1]] , [np.log10(np.log10(m))]*2, '-', c=c[ii], lw=1.5, alpha=0.9 ) 
            
            ax1.text(ax1.get_xlim()[0]+0.3, np.log10(np.log10(m))+0.02, 'Mach$=$%4.1f' % (m), fontsize=10, color=c[ii], alpha=0.9)
            ax2.text(ax2.get_xlim()[0]+0.3, np.log10(np.log10(m))+0.02, 'Mach$=$%4.1f' % (m), fontsize=10, color=c[ii], alpha=0.9)
            
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im2, cax=cbar_ax)
        plt.savefig(saveplot + '.png')
        
    if plot_new:
        fig = plt.figure()
        
        ax1 = plt.subplot(111)
        im2 = ax1.pcolor( np.log10(keV2K) + xnew, ynew, Znew, vmin=-8)
#        ax1.set_title("Interpolated function")
        ax1.set_xlim([7, 8.4])
        ax1.set_ylim([np.log10(np.log10(1.7)), np.log10(np.log10(10.))])
        ax1.set_xlabel('$\\mathrm{log_{10}(T)\\,[K]}$ ')
        ax1.set_ylabel('$M$ ')
        
        y_ticks = [np.log10(np.log10(m)) for m in [1.7,2.5,4,10]]
        print ['%.2e' % (y) for y  in y_ticks], [10**(10**y) for y in y_ticks]
        ax1.set_yticklabels([10**(10**y) for y in y_ticks])
        plt.yticks(y_ticks)
    
#        temp  = [1.5,2.2,3.0,10.0]
#        c     = [plt.cm.rainbow( (np.log10(np.log10(m))-ax1.get_ylim()[0])/abs(ax1.get_ylim()[1]-ax1.get_ylim()[0]) ) for m in mach]
#        for ii,m in enumerate(mach):
#            ax1.plot( [ax1.get_xlim()[0], ax1.get_xlim()[1]] , [np.log10(np.log10(m))]*2, '-', c=c[ii], lw=1.5, alpha=0.9 ) 
#            ax2.plot( [ax2.get_xlim()[0], ax2.get_xlim()[1]] , [np.log10(np.log10(m))]*2, '-', c=c[ii], lw=1.5, alpha=0.9 ) 
#            
#            ax1.text(ax1.get_xlim()[0]+0.3, np.log10(np.log10(m))+0.02, 'Mach$=$%4.1f' % (m), fontsize=10, color=c[ii], alpha=0.9)
#            ax2.text(ax2.get_xlim()[0]+0.3, np.log10(np.log10(m))+0.02, 'Mach$=$%4.1f' % (m), fontsize=10, color=c[ii], alpha=0.9)
            
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im2, cax=cbar_ax, label='$\log_{10}\Phi$',)
        plt.savefig(saveplot + '_DSA.pdf')
        plt.savefig(saveplot + '_DSA.png', dpi=800)       
        
        
    if plot_PhD:
        fig = plt.figure()
        
        temp = np.linspace(2,20,20)
        print temp
        mach = np.linspace(2,7,300)
        psi_x,psi_y = LoadFile_psi(psiFile)
        import itertools
        H,M,T  = [],[],[]
        for t in temp:
            results_temp = math.find_closest(psi_x, t)
            results_mach = math.find_closest(psi_y, mach) # 
            H.append(H_psi[results_mach,np.ones_like(results_mach)*results_temp])
            M.append(mach)
            T.append(np.ones_like(results_mach)*t)
           
        H = list(itertools.chain.from_iterable(H))
        M = list(itertools.chain.from_iterable(M))
        T = list(itertools.chain.from_iterable(T))  
        
        plt.scatter(M,np.log10(H),c=T,alpha=0.1,s=5) 
        cb = plt.colorbar(label='Downstream Temperature [keV]')
        cb.set_alpha(1)
        cb.draw_all()
        plt.xlabel('Mach number $M$')
        plt.ylabel('$\log_{10}\,\Phi(M,T)$')
        plt.savefig(saveplot + '_PhD.pdf')
        plt.savefig(saveplot + '_PhD.png', dpi=800)       
        
        
    # Save File A
    if psiFileNew:
      location = psiFileNew
    else:
      location = psiFile.replace('.txt', 'fine(%i,%i).txt' % (inter[0],inter[1]) )  
      
    header = '#    Mach'
    for x in xnew:
        header += '%13.4e' % (10**x)

    mf        = open(location,"w")
    mf.write(header + '\n')
   
    for ii,y in enumerate(ynew):           
        string = '%9.4f' % (10**(10**y)) +  ''.join(['%13.4e' % (10**z) for z in Znew[ii][:]])
        mf.write(string + '\n') 

    mf.close()  
    
    
    #=== FILE B ===#
    Y_new    = np.empty( (1,1) )
    for ii,h in  enumerate(H_mach.T):
       interp_spline = interpolate.interp1d( 10**psi_y , h, kind='cubic') 
       #print np.expand_dims(interp_spline( 10**ynew ), axis=1).shape
       
       if Y_new.shape[0] > 1:
          Y_new         = np.hstack( (Y_new, np.expand_dims(interp_spline( 10**ynew ), axis=1) )   )
       else:
          Y_new         = np.expand_dims(interp_spline( 10**ynew ), axis=1)
       #print Y_new.shape
       
    # Save File B
    if machFileNew:
      location = machFileNew
    else:
      location = machFile.replace('.txt', 'fine(%i,%i).txt' % (inter[0],inter[1]) )  
      
    header = '#         q            M           r           M*(1-1/r)        s'
     
    mf        = open(location,"w")
    mf.write(header + '\n')

    for ii,y in enumerate(10**ynew):    
        string = ''.join(['%14.6e' % (y) for y in Y_new[:][ii]]) #some numbers are very large and ewould need a good margin
        mf.write(string + '\n') 

    mf.close()  

    return 0
    

if __name__ == "__main__":
    InterpolateRadio2D(psiFile = '../Analysis_MUSIC2/Hoeft_radio/mach_psi_tablefine(10,3).txt',  inter=(3,2)) #(90,27)  

