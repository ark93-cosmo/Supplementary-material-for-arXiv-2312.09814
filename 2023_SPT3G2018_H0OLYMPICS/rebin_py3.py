import numpy as nm
import copy
def get_l_rebinning(lcuts, lmin, lmax,ellb=None, Dl=False):#rebinning must be subsample of the original binning, ellb are the binned ell's
#    if lmin >= lmax:
#        bins = None
#        return bins
    if isinstance(lcuts, str):
        lcuts = nm.loadtxt(lcuts)
#    if lcuts is None:
#        lcuts = arange(lmin, lmax+1)
    imin = lcuts.searchsorted(lmin)
    if lcuts[imin] != lmin:
        lcuts = nm.hstack(([lmin], lcuts[imin:]))
    else:
        lcuts = lcuts[imin:]
    imax  = lcuts.searchsorted(lmax+1)
    lcuts = nm.hstack((lcuts[:imax], [lmax+1]))
    nbins = len(lcuts) - 1
#    bins = zeros([nbins, lmax + 1 - lmin])
    bins = nm.zeros([nbins, len(ellb)])
   

        
    for i, lp in enumerate(zip(lcuts[:-1], lcuts[1:])):
        lm  = lp[0] - lmin 
        lM  = lp[1] - 1 - lmin

        if any(ellb)==None:
            ell = arange(lm+lmin, lM + 1 + lmin)    
            bins[i,lm:lM+1] = ell * (ell + 1.0)
            nrm = sum(bins[i,lm:lM+1])
            bins[i,lm:lM+1]/=nrm
        else:

            index=nm.array([ii for ii in range(len(ellb))if  ellb[ii]>=lm+lmin and ellb[ii]< lM + 1 + lmin])
            ell=copy.deepcopy(ellb[index])
            if Dl:
                bins[i,index[0]:index[-1]+1] = ell * (ell + 1.0)*0.+1. #Do not give weight if I bin Dls                
            else:
                bins[i,index[0]:index[-1]+1] = ell * (ell + 1.0)
            nrm = sum(bins[i,index[0]:index[-1]+1])
            bins[i,index[0]:index[-1]+1]/=nrm
    
    return bins

def rebin_cov(cov, bin_total):
    if bin_total.shape[0] == bin_total.shape[1]:
        return cov
    cov = nm.dot(bin_total, cov)
    cov = nm.dot(cov, bin_total.T)
    return cov

def rebin_vec(vec, bin_total):
    if bin_total.shape[0] == bin_total.shape[1]:
        return vec
    vec = nm.dot(bin_total, vec)

    return vec


def rebin_diag_cov(err,bin_total,var=True):
    cov=nm.identity(len(err))
#    print len(err)                                                                                                                       
    if not var:
        err=err**2
#        print 'not true'                                                                                                                 
    for ii in range(len(err)):
        cov[ii,ii]=err[ii]
    where_are_NaNs = nm.isnan(cov)
    cov[where_are_NaNs]=0
    return rebin_cov(cov, bin_total)


def plot_ps_resid_split(title,lm0,yps0,errps0,modelps0,yres0,errres0,lmin,lmax,bsplit,Cl=False,lmunbin=None,modelpsunbin=None,reslims=None):
    #plt.ion()
    
    cnavy='#000080'
    fonts=18
    fontstick=fonts-3
    lw=1.5
    fig=plt.figure(figsize=(8,10))
    f, (ax1,ax2) = plt.subplots(2, sharex=True, figsize=(8,10))
    f.subplots_adjust(left=0.14, right=0.9, bottom=0.09, top=0.93,hspace=0)
    f.suptitle(title, y=0.96, fontsize=fonts-5)

    #First plot:power spectrum left
    lm=copy.deepcopy(lm0)
    factlm=lm0*(lm0+1)/2/nm.pi if not Cl else 1.
    modelps=copy.deepcopy(modelps0)
    
    if((modelpsunbin==None and not lmunbin==None) or (not modelpsunbin==None and lmunbin==None  )):
        print('you must provide both unbinned lm and modelps for your model')
        stop()
    if (modelpsunbin==None and  lmunbin==None): 
        ax1.plot(lm,modelps*factlm,'r-')
    else:
        factlmunbin=lmunbin*(lmunbin+1)/2/nm.pi if not Cl else 1.
        ax1.plot(lmunbin,modelpsunbin*factlmunbin,'r-')
    limsmod=ax1.axis()

    temp=yps0*factlm
    tempmin=(modelps0-errps0)*factlm
    tempmax=(modelps0+errps0)*factlm
    
    bmin=nm.argmin(tempmin[:bsplit])
    bmax=nm.argmax(tempmax[:bsplit])
    bminR=nm.argmin(temp[bsplit:])  #if bsplit<len(lm0)-1 else 0
    bmaxR=nm.argmax(temp[bsplit:]) #if bsplit<len(lm0)-1 else 0
    #fact=max(abs(temp[bsplit+bminR]/tempmin[bmin]),temp[bsplit+bmaxR]/tempmax[bmax])

    yps=copy.deepcopy(yps0)
    errps=copy.deepcopy(errps0)

    factlm=lm*(lm+1)/2/nm.pi if not Cl else 1.

    yy=copy.deepcopy(modelps)
    c1d='#C0C0C0'
    c2d='#DCDCDC'
    c3d='#f0f0f0'
    #ax1.plot(lm,(yy+errps)*factlm,color=c1d)
    #ax1.plot(lm,(yy-errps)*factlm,color=c1d)
    #ax1.plot(lm,(yy+2*errps)*factlm,color=c2d)
    #ax1.plot(lm,(yy-2*errps)*factlm,color=c2d)

    #ax1.fill_between(lm,(yy-3*errps)*factlm, (yy+3*errps)*factlm,color=c3d)
    #ax1.fill_between(lm,(yy-2*errps)*factlm, (yy+2*errps)*factlm,color=c2d)
    #ax1.fill_between(lm,(yy-errps)*factlm, (yy+errps)*factlm,color=c1d)

    #ax1.plot(lm,yps*factlm,color=cnavy,marker='o',ms=2,linestyle='-')
    ax1.errorbar(lm,yps*factlm,yerr=errps*factlm,color=cnavy,marker='o',ms=2,linestyle='None',capthick=0)
    
    ax1.set_xlabel('$\ell$', fontsize=fonts)
    if Cl:
        ax1.set_ylabel('$ \mathcal{C}_{\ell} \ [\mathrm{\mu K}^2]$', fontsize=fonts)
    else:
        ax1.set_ylabel('$ \mathcal{D}_{\ell} \ [\mathrm{\mu K}^2]$', fontsize=fonts)
    

    ax1.tick_params(axis='both', which='major', labelsize=fontstick)

    ax1.set_ylim(min(limsmod[2],tempmin[bmin]),max(limsmod[3],tempmax[bmax]*1.1))
    
    ax1.minorticks_on()
    ax1.set_xlim(lmin,lmax)
 
    #Second plot:residuals left
    ax2.plot(lm,lm*0.,'r')

    factlm=lm0*(lm0+1)/2/nm.pi if not Cl else 1.
    temp=yres0*factlm
    bmin=nm.argmin((temp)[:bsplit])
    bmax=nm.argmax((temp)[:bsplit])
    bminR=nm.argmin((temp)[bsplit:]) #if bsplit<len(lm0)-1 else 0
    bmaxR=nm.argmax((temp)[bsplit:]) #if bsplit<len(lm0)-1 else 0

    #fact=max(abs(temp[bsplit+bminR]/temp[bmin]),temp[bsplit+bmaxR]/temp[bmax])

    fact=max(abs(temp[bsplit+bminR]/temp[bmin]),temp[bsplit+bmaxR]/temp[bmax])
    
    lm=lm0[:bsplit]
    yres=yres0[:bsplit]
    errres=errres0[:bsplit]
    factlm=lm*(lm+1)/2/nm.pi if not Cl else 1.
    yy=0 #yres #0
    #ax2.plot(lm,yy+errres*factlm,color=c1d)
    #ax2.plot(lm,yy-errres*factlm,color=c1d)
    #ax2.plot(lm,(yy+2*errres)*factlm,color=c2d)
    #ax2.plot(lm,-(yy+2*errres)*factlm,color=c2d)

#    ax2.fill_between(lm,yy-3*errres*factlm, yy+3*errres*factlm,color=c3d)
#    ax2.fill_between(lm,yy-2*errres*factlm, yy+2*errres*factlm,color=c2d)
#    ax2.fill_between(lm,yy-errres*factlm, yy+errres*factlm,color=c1d)

 #   ax2.plot(lm,yres*factlm,color=cnavy,marker='o',ms=2,linestyle='-')
    ax2.errorbar(lm,yres*factlm,yerr=errres*factlm,color=cnavy,marker='o',ms=2,linestyle='None',capthick=0)
 

    # Second plot:residuals rigth


    ax2.plot(lm,lm*0.,'r')
    if bsplit==len(lm0)-1: fact=1
    lm=lm0[bsplit:]
    yres=yres0[bsplit:]/fact
    errres=errres0[bsplit:]/fact
    factlm=lm*(lm+1)/2/nm.pi if not Cl else 1.
    yy=0 #yres
    # ax2.plot(lm,yy+errres*factlm,color=c1d)
    # ax2.plot(lm,yy-errres*factlm,color=c1d)
    # ax2.plot(lm,(yy+2*errres)*factlm,color=c2d)
    # ax2.plot(lm,-(yy+2*errres)*factlm,color=c2d)
    #    ax2.fill_between(lm,yy-3*errres*factlm, yy+3*errres*factlm,color=c3d)
    #    ax2.fill_between(lm,yy-2*errres*factlm, yy+2*errres*factlm,color=c2d)
    #    ax2.fill_between(lm,yy-errres*factlm, yy+errres*factlm,color=c1d)
    
    #    ax2.plot(lm,yres*factlm,color=cnavy,marker='o',ms=2,linestyle='-')
    ax2.errorbar(lm,yres*factlm,yerr=errres*factlm,color=cnavy,marker='o',ms=2,linestyle='None',capthick=0)

    
    factlm=lm0[bmin]*(lm0[bmin]+1)/2/nm.pi if not Cl else 1.
    factlm2=lm0[bmax]*(lm0[bmax]+1)/2/nm.pi if not Cl else 1.
    ax2lim=max(abs(yres0[bmin]*factlm),abs(yres0[bmax]*factlm2))
    # ax2.set_ylim(yres0[bmin]*factlm,yres0[bmax]*factlm2)
    
    ax2.set_ylim(-ax2lim,ax2lim) if reslims==None else ax2.set_ylim(reslims[0],reslims[1])    
    if not bsplit==len(lm0)-1:    
        ax2.plot([lm0[bsplit],lm0[bsplit]],[-1.e12,1.e12],':k',lw=lw)    
        ay2 = ax2.twinx()
        factlm2=lm0[bmax]*(lm0[bmax]+1)/2/nm.pi if not Cl else 1.
        ay2lim=max(abs(yres0[bmin]*fact*factlm),abs(yres0[bmax]*fact*factlm2))
        # ay2.set_ylim(yres0[bmin]*fact*factlm,yres0[bmax]*fact*factlm2)
        ay2.set_ylim(-ay2lim,ay2lim)
        ay2.tick_params(axis='y', which='major', labelsize=fontstick)

    
    ax2.minorticks_on()
    ax2.set_xlabel('$\ell$', fontsize=fonts)
    if Cl:
        ax2.set_ylabel('$\Delta \mathcal{C}_{\ell} \ [\mathrm{\mu K}^2]$', fontsize=fonts)
    else:
        ax2.set_ylabel('$\Delta \mathcal{D}_{\ell} \ [\mathrm{\mu K}^2]$', fontsize=fonts)
    ax2.tick_params(axis='both', which='major', labelsize=fontstick)
    ax2.set_xlim(lmin,lmax)
    plt.tight_layout()
