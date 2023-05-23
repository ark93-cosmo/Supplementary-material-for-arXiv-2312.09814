
import numpy as nm
try:
    import pyfits as pf
except ImportError:
    import astropy.io.fits as pf
import os
import copy
from rebin_py3 import *

def mkdirsafe(dirs):
    import os
    if not os.path.exists(dirs):  os.makedirs(dirs)
    return


def nonans(invec):
    wherenan=nm.isnan(invec)
    invec[wherenan]=0
    return invec

def tick_function(X):
    tot=[]
    theta=180./nm.array(X)
    for z in theta:
        if z>=1:
            tot.append('%.f' % z+'$^\circ$')
        elif z>=0.1:
            tot.append('%.1f' % z+'$^\circ$')
        elif z>=0.01:
            tot.append('%.2f' % z+'$^\circ$')

    return tot

def cm2inch(cm):
    """Centimeters to inches"""
    return cm *0.393701


class class_spectra_properties:
    def __init__(self,Cl):
        self.Cl=Cl

class class_plot_bleaks:
    def __init__(self,plot_residuals):

        if plot_residuals:
            self.plot_bleaks=[ False] #[True,False]
        else:
            self.plot_bleaks=[False]

class class_colors:
    def __init__(self):
        self.colord='gray'
        #classical
        cols={'cnavyl':'#ff4709','cnavy':'#000080','cnavyCV':'#a6d2ff','colorbf':'#0984ff'}
        #roma=
        #            cols={'cnavy':'k','cnavyl':'#a7051b','cnavyCV':'#ffd269','colorbf':'#F2A900'}
        #europe
        #            cols={'cnavy':'k','cnavyl':'#001489','cnavyCV':'#ffd269','colorbf':'#FEDD00'}
        
        self.cnavy=cols['cnavy']# 'k' #'#000080'
        self.cnavys=[cols['cnavy'],'#d3d3d3']# 'k' #'#000080'
        
        #            cnavyl='#00a180' # '#0984ff'
        #            cnavyCV='#c8fff4' #'#a6d2ff'
        
        self.cnavyl= cols['cnavyl'] #'#F2A900' #'#ff4709' #'r' # '#bf00ef' #'#30a9ff' #'#3098ff' #'#4a4aff'#,'#3030ff'
        self.cnavyls= [cols['cnavyl'],'#6495ed']
        self.cnavyCV=cols['cnavyCV'] #'#a6d2ff'
        self.colbin1err='#FFA500'
        self.colbin1='#c1c1c1'
        self.colorbf=cols['colorbf'] #'#0984ff'
        #Giallo Roma (130C)#F2A900 e il Rosso Roma (202C)#862633 , cv #ffd269
        #european flag 001489 pantone reflex blue C, pantone yellow FEDD00
        

class class_plot_properties:
    def __init__(self,plot_bleak,width,plot_residuals,plot_unbinned,polstoplot0):
        
        self.lwd=0.5
        self.dashes = [2,2,2,2]
        self.labpad=10
        self.labpadCL=10
        
        if plot_bleak :
            self.polstoplot=['EEd','TE','EE']
        else:
            self.polstoplot=polstoplot0

        if  plot_residuals or plot_unbinned:
            self.plot_cvl=False
        else:
            self.plot_cvl=True

        if not plot_residuals:
            self.topaxis=False #True
        else:
            self.topaxis=False
        if width==18:
            self.lw=1.1
            self.lw1=0.8
            self.labelx0=-0.08
            self.labelxCL=-0.1
        else:
            self.lw=0.8
            self.lw1=0.8
            self.labelx0=-0.1
            self.labelxCL=-0.13


        self.lt=5 if width==18 else 3
        self.ltm=3 if width==18 else 2
        self.lwt=0.8 if width==18 else 0.5 #width of the ticks
        if self.topaxis:
            self.lwt=1. if width==18 else 0.8
        if plot_residuals:
            self.fonts=13 if width==18 else 9
            self.fonts1=11 if width==18 else 8
        else:
            self.fonts=14 if width==18 else 8
            self.fonts1=12 if width==18 else 7

        self.mss0=3 if width==18 else 2.2
        self.mssunbinned=0.8 if width==18 else 0.6
         
#def change_plot_properties_pol(pol,ppr):
    

    
class class_ticks():
    def __init__(self,plot_unbinned,leftpols,highl=[],lowl={}):
        self.Clfact={}
        self.yticks={}
        self.yticksres={}
        self.ytickslab={}
        self.yticksreslab={}
        self.yticksresleft={}
        self.yticksleftPS={}
        self.yticksleftPSlab={}
        self.ylimsleftPS={}           
        
        self.logscales=['PP']
        # x limits and ticks
        self.xlims={'TT':[1.5,2550],'EE':[10,2050],'TE':[10,2050],'PP':[8,2048]}
        self.xticksR={'TT':[30,500,1000,1500,2000,2500],'TE':[30,500,1000,1500,2000],'EE':[30,500,1000,1500,2000],'PP':[10,100,1000]}
        self.xticksL={'TT':[2,10],'TE':[2,10],'EE':[2,10]}

        #ylimits
        self.ylims={'TT':[-200,6200],'EEd':[-10,45],'TE':[-170,170],'EEc':[-0.0001, 0.0011],'PP':[0,1.8e-7]}
        self.ylimsres={'TT':[-80,80],'EEd':[-2,2],'TE':[-15,15],'EEc':[-0.00006, 0.00006],'PP':[-0.5e-7,0.5e-7]}
        self.ylimsresleft={'TT':[-700,700],'TE':[-20,20],'EEd':[-0.7,0.7],'EEc':[-0.001, 0.001]}


        #yticks
        self.yticks['TT']=[0,1000,2000,3000,4000,5000,6000]
        self.yticksresleft['TT']=[-600,-300,0,300,600]
        self.yticksresleft['EEd']=[-0.6,-0.3,0,0.3,0.6]
        self.yticksresleft['EEc']=[-0.05,0,0.05]
        self.yticksresleft['TE']=[-16,-8,0,8,16]

        
        self.yticksres['TT']=[-60,-30,0,30,60] 
        self.yticks['EEd']= [0,10,20,30,40]
        self.yticksres['EEd']=[-1,0,1]
        self.yticks['TE']=[-140,-70,0,70,140]
        self.yticksres['TE']=[-10,0,10]

        self.yticks['PP']= [2e-8,4e-8,6e-8,8e-8,10e-8,12e-8,14e-8,16e-8]
        self.yticksres['PP']=[-0.4e-7,-0.2e-7,0,0.2e-7,0.4e-7]

        #ylimsleftPS['EEd']=[-0.3,0.3]
        #yticksleftPS['EEd']=[aa*1.e-4 for aa in [0,0.3,0.1]]

        #Limits for EE Cl
#        self.ylims[#0.001  #[aa*1.e-4 for aa in range(0,10)]#[-0.0001,0.001]

#        self.ylimsres
        self.yticks['EEc']=[aa*1.e-4 for aa in range(0,12,2)] #[ylims['EE'][0],0,ylims['EE'][1]]
        self.yticksres['EEc']=[-0.00004,0,0.00004]


        if plot_unbinned:
            self.yticksres['TT']=[-200,-100,0,100,200] 
            self.ylimsres['TT']=[-250,250]
            self.yticksres['TE']=[-25,0,25]
            self.ylimsres['TE']=[-30,30]
            self.yticksres['EEd']=[-3,0,3]
            self.ylimsres['EEd']=[-4,4]
            self.yticksres['EEc']=[-0.00005,0,0.00005]
            self.ylimsres['EEc']=[-0.0001, 0.0001]

        #            ylimsleftPS['EEc']=ylims['EEc']
        #            yticksleftPS['EEc']= yticks['EEc']
        self.ylimsleftPS['EEc']=[-0.05,0.2]
        self.yticksleftPS['EEc']=[aa for aa in [0,0.05,0.1,0.150005]]

        self.ylimsleftPS['EEd']=[-0.65,3]
        self.yticksleftPS['EEd']=[0,1,2]

        self.ylimsleftPS['TE']=[-25,28]
        self.yticksleftPS['TE']=[-20,-10,0,10,20]



        if 'EE' in lowl.keys() or 'TE' in lowl.keys():
            if 'HFI' in lowl['EE']:

                self.ylimsleftPS['EEc']=[-0.05,0.2]
                self.yticksleftPS['EEc']=[aa for aa in [0,0.05,0.1,0.150005]]
                
                self.ylimsleftPS['EEd']=[-0.03,0.23]
                self.yticksleftPS['EEd']=[0,0.1,0.2]

                
                self.ylimsleftPS['TE']=[-6.3,10]
                self.yticksleftPS['TE']=[-4,0,4,8]

                self.ylimsresleft['TE']=[-10,10]
                self.yticksresleft['TE']=[-8,-4,0,4,8]
                                
                self.ylimsresleft['EEd']=[-0.16,0.16]
                self.yticksresleft['EEd']=[-0.1,0,0.1]
                
                self.ylimsresleft['EEc']=[-0.15, 0.15]
                self.yticksresleft['EEc']=[-0.1,0,0.1]

                

        #minorticks
        self.minorrangelog=range(2,10,1)+range(10,30,10) #[2,3,4,5,6,7,8,9,20]
        self.minorrangelin={}

        for pol in ['TT','TE','EE']:
            self.minorrangelin[pol]=range(100,self.xlims[pol][1],100) #[2,3,4,5,6,7,8,9,20]
        self.minorrangelin['PP']=range(self.xlims['PP'][0],10,1)+range(10,100,10)+range(100,1000,100)+range(1000,self.xlims[pol][1],1000) #[2,3,4,5,6,7,8,9,20]
        def findell(ell):
            return 180./nm.array(ell)
            
        self.xticksRTOP={}
        self.xticksRTOP['TT']=findell([1,0.2,0.1,0.07])
        self.xticksRTOP['EE']=findell([1,0.2,0.1])
        self.xticksRTOP['TE']=findell([1,0.2,0.1])



    
        ###################################

        #Set variables for lin-log plot
        ###################################

        self.xminlog=self.xlims['TT'][0]
        self.xmid=30.

        self.xminlin={}
        self.xminlin['TT']=-800. if 'TT' in leftpols else self.xlims['TT'][0]
        self.xminlin['EE']=-600. if 'EE' in leftpols else self.xlims['EE'][0]
        self.xminlin['TE']=-600. if 'TE' in leftpols else self.xlims['TE'][0]
        self.xminlin['PP']=self.xlims['PP'][0]
        
        self.xmaxlin={}
        self.xmaxlin['TT']=self.xlims['TT'][1]
        self.xmaxlin['TE']=self.xlims['TE'][1]
        self.xmaxlin['EE']=self.xlims['EE'][1]
        self.xmaxlin['EE']=self.xlims['PP'][1]

def change_tks(pol,Cl,tks,ALLlow):
    tks.Clfact={'EEc':1e-5}
    tks.Clfact['PP']=1e-8
    if pol=='EEd':
        Cl['EE']=False
        pol='EE'
    else:
        Cl['EE']=True

    if Cl['EE']:
        tks.ylims['EE']=tks.ylims['EEc']
        tks.ylimsres['EE']=tks.ylimsres['EEc']
        
        tks.Clfact['EE']=tks.Clfact['EEc']
        tks.yticks['EE']=tks.yticks['EEc']
        tks.yticksres['EE']=tks.yticksres['EEc']
        tks.ylimsresleft['EE']=tks.ylimsresleft['EEc']
        tks.yticksresleft['EE']=tks.yticksresleft['EEc']
        tks.ylimsleftPS['EE']=tks.ylimsleftPS['EEc']
        tks.yticksleftPS['EE']=tks.yticksleftPS['EEc']

    else:
        tks.ylims['EE']=tks.ylims['EEd']
        tks.ylimsres['EE']=tks.ylimsres['EEd']
#        tks.Clfact={}
        tks.yticks['EE']=tks.yticks['EEd']
        tks.yticksres['EE']=tks.yticksres['EEd']
        tks.ylimsresleft['EE']=tks.ylimsresleft['EEd']
        tks.yticksresleft['EE']=tks.yticksresleft['EEd']
        tks.ylimsleftPS.pop("EE", None)
        if ALLlow:
            tks.ylimsleftPS['EE']=tks.ylimsleftPS['EEd']
            tks.yticksleftPS['EE']=tks.yticksleftPS['EEd']
    return tks,pol




class class_experiment():

    def __init__(self,filename,errfilename, modelfile,calplanck):
        self.binningfile='binall.dat'
        self.relcuts='binall.dat'
        self.filename=filename
        self.errfilename=errfilename
        self.model=nm.loadtxt(modelfile)
        #To check: now all the models for all data are rescaled by calplanck, and not just plik
        for ii in range(1,6):self.model[:,ii]=self.model[:,ii]/calplanck**2
        self.lm={}
        self.lm['TT']=[2,29]
        self.lm['TE']=[2,29]
        self.lm['EE']=[2,29]
        self.indmod={'TT':1,'TE':2,'EE':3}
        self.fsky={}
     #   print 'shape mode', self.model[:,ii].shape

def read_covmat(errfilename):
        print ('read covmat')
        # read covmat
        #print 'load covmat'
        # Read binary if it exists, faster
        
        if os.path.exists((errfilename).replace('.dat','.binary')):
            covmat=nm.fromfile((errfilename).replace('.dat','.binary'))
        else:
            covmat=nm.loadtxt(errfilename)
        f=open(errfilename)
        #read shape of covmat in first line
        covshape=[float(aa) for aa in f.readline().split('shape: (')[-1].split(')')[0].split(',')]
        f.close()
        covmat.shape=(covshape[0],covshape[1])
        #print 'Covmat read plik'
        return covmat

class class_plik(class_experiment):

    def __init__(self,filename,errfilename, modelfile,calplanck,is_bin1=True):
        class_experiment.__init__(self,filename,errfilename, modelfile,calplanck)
        if is_bin1:
            self.binningfile='binall.dat' #unbinned
            self.relcuts='bin30.dat'      #delta=30
        else:
            self.binningfile='binning.dat' #Plik binning
            self.relcuts='rebinning30.dat' #Approx rebinning delta=30 of the Plik binning
        #read power
        self.power=nm.loadtxt(self.filename)
        #indexes

        self.ind={'TT':1,'EE':3,'TE':5} #,'errpTT':2,'errnTT':2,'errpTE':2,'errnTE':2,'errpEE':2,'errnEE':2}
        #limits and fsky
        self.lm['TT']=[30,2508]
        self.lm['TE']=[30,1996]
        self.lm['EE']=[30,1996]

        self.fsky['TT']=0.66
        self.fsky['EE']=0.70
        self.fsky['TE']=0.70
        self.is_bin1=is_bin1

    def def_polb(self):
        polb={}
        allcov=True #If you have the full TT,TE, EE covmat
        if allcov:
            if self.is_bin1:
                polb['TT']=nm.array([0,2508-30+1])
                polb['EE']=polb['TT'][-1]+nm.array([0,1996-30+1])
                polb['TE']=polb['EE'][-1]+nm.array([0,1996-30+1])
            else:
                polb['TT']=nm.array([0,215])
                polb['EE']=polb['TT'][-1]+nm.array([0,199])
                polb['TE']=polb['EE'][-1]+nm.array([0,199])
        else:
            if self.is_bin1:
                polb['TT']=nm.array([0,2508-30+1])
                polb['EE']=nm.array([0,1996-30+1])
                polb['TE']=nm.array([0,1996-30+1])
        return polb
    def read_and_bin_covmat(self,pols,name=None):
        print ('read covmat and bin')
        polb=def_polb(self)
                #Define indexes max-min to cut the covmat in polarizations
        covmat=read_covmat(self.errfilename)
        self.covmat=covmat
        polb=def_polb(self)
        errbpp={}
        errbpn={}

        errbpp['plik']={}
        
        for pol in pols:
            #print 'start binning', pol
            # If the starting comvat is binned, I first reconstruct the correspondin ll for each frequency and then find the correct bins by rebinning this.
            lltemp=nm.array(nm.arange(3000))
            if not self.binningfile=='binall.dat':
                binstemp=get_l_rebinning(self.binningfile,self.lm[pol][0],self.lm[pol][1],ellb=lltemp)
                #print 'binning matrix '
                lmcov=rebin_vec(lltemp,binstemp)
            else:
                lmcov=lltemp[self.lm[pol][0]:self.lm[pol][1]+1]
            covpol=covmat[polb[pol][0]:polb[pol][1],:][:,polb[pol][0]:polb[pol][1]]
            #print 'polb===='
        
            #print polb[pol][0],polb[pol][1]
            binstemp=get_l_rebinning(self.relcuts,self.lm[pol][0],self.lm[pol][1],ellb=lmcov,Dl=False)

            blmcov=rebin_vec(lmcov,binstemp)
            #print 'start covmat binning'
            errbpp['plik'+pol]=nm.sqrt(nm.diag(rebin_cov(covpol,binstemp)))*blmcov*(blmcov+1)/2./nm.pi
            #print 'done'
            errbpn['plik'+pol]=copy.deepcopy(errbpp['plik'+pol])
            # Memorize as well the unbinned diagonals
            errbpp['plik'+pol+'bin1']=nm.sqrt(nm.diag(covpol))*lmcov*(lmcov+1)/2./nm.pi
            errbpn['plik'+pol+'bin1']=copy.deepcopy(errbpp['plik'+pol+'bin1'])
            errbpp['plik'+pol+'bin1'+'lm']=lmcov
        return errbpp,errbpn



            
class class_commander(class_experiment):
    def __init__(self,filename,errfilename, modelfile,calplanck):
        
        class_experiment.__init__(self,filename,errfilename, modelfile,calplanck)
        self.power=nm.loadtxt(self.filename)

        self.ind={'TT':3,'errpTT':4,'errnTT':5}
        self.fsky['TT']=0.93



class class_LFI(class_experiment):
    def __init__(self,filename,errfilename, modelfile,calplanck):

        class_experiment.__init__(self,filename,errfilename, modelfile,calplanck)

        temp=pf.open(self.filename)[1].data.astype([('TEMPERATURE', '<f8'), ('GRADIENT', '<f8'), ('CURL', '<f8'), ('G-T', '<f8'), ('C-T', '<f8'), ('C-G', '<f8')]).view('<f8')
        err=pf.open(self.errfilename)[1].data.astype([('TEMPERATURE', '<f8'), ('GRADIENT', '<f8'), ('CURL', '<f8'), ('G-T', '<f8'), ('C-T', '<f8'), ('C-G', '<f8')]).view('<f8')
        temp.shape=(30,6)
        lllfi=nm.arange(30)
        lllfi.shape=(30,1)
        err.shape=(30,6)
        self.power=nm.hstack((lllfi,temp,err))
                                                                                              

        self.ind={'TT':1+0,'EE':1+1,'TE':1+3,'errpTT':1+0+6,'errpEE':1+1+6,'errpTE':1+3+6,'errnTT':1+0+6,'errnEE':1+1+6,'errnTE':1+3+6}
        self.fsky['EE']=1.
        self.fsky['TE']=1.
#        print lllfi
        print ('shape lFI', self.power.shape)
class class_HFI(class_experiment):

    def __init__(self,filename,errfilename, modelfile,calplanck):
        class_experiment.__init__(self,filename,errfilename, modelfile,calplanck)
        lm= nm.loadtxt((self.filename).replace('XX','EE'))[:,0]
        fact=1. #lm*(lm+1.)/nm.pi/2.
        #Luca's files are Dl. power should be in Dl
        tempEE=nm.loadtxt(self.filename.replace('XX','EE'))[:,1]*fact
        errtempEE=nm.loadtxt(self.filename.replace('XX','EE'))[:,2]*fact

        tempTE=nm.loadtxt(self.filename.replace('XX','TE'))[:,1]*fact
        errtempTE=nm.loadtxt(self.filename.replace('XX','TE'))[:,2]*fact

        lm=nm.hstack((0,1,lm))
        
        tempEE=nm.hstack((0,0,tempEE))
        errtempEE=nm.hstack((0,0,errtempEE))
        tempTE=nm.hstack((0,0,tempTE))
        errtempTE=nm.hstack((0,0,errtempTE))
        #in DL!
        self.power=nm.transpose(nm.vstack((lm,tempEE,errtempEE,tempTE,errtempTE)))
        print ('shape', self.power.shape)
        
        #        self.ind={'TT':1+0,'EE':1+1,'TE':1+3,'errpTT':1+0+6,'errpEE':1+1+6,'errpTE':1+3+6,'errnTT':1+0+6,'errnEE':1+1+6,'errnTE':1+3+6}


        self.ind={'EE':1,'errpEE':2,'errnEE':2,'TE':3,'errpTE':4,'errnTE':4}
        self.fsky['EE']=1

class class_lensing():

    def __init__(self,filename,errfilename, modelfile,calplanck,binnedmodel):
        #        class_experiment.__init__(self,filename,errfilename, modelfile,calplanck)
        self.filename=filename
        self.binningfile=nm.loadtxt((self.filename))[:,1]
        self.relcuts=nm.loadtxt((self.filename))[:,1]

        self.errfilename=errfilename
        self.model=nm.loadtxt(modelfile)
        #To check: now all the models for all data are rescaled by calplanck, and not just plik
        for ii in range(1,6):self.model[:,ii]=self.model[:,ii] #TBCHECKED
        self.lm={}
        self.lm['PP']=[8,2048]
        self.indmod={'TT':1,'TE':2,'EE':3,'PP':5}
        self.fsky={}

        lm= nm.loadtxt((self.filename))[:,3]
        fact=lm*(lm+1.)/nm.pi/2.
        temp=nm.loadtxt(self.filename)[:,4]
        err=nm.loadtxt(self.filename)[:,5]

        self.power=nm.transpose(nm.vstack((lm,temp,err)))
        print ('shape', self.power.shape)
        
        #        self.ind={'TT':1+0,'EE':1+1,'TE':1+3,'errpTT':1+0+6,'errpEE':1+1+6,'errpTE':1+3+6,'errnTT':1+0+6,'errnEE':1+1+6,'errnTE':1+3+6}


        self.ind={'PP':1,'errpPP':2,'errnPP':2}
        self.fsky['PP']=1
        self.binnedmodel=nm.loadtxt(binnedmodel)

class class_lensing2():

    def __init__(self,filename,errfilename, modelfile,calplanck,binnedmodel):
        #this is not good since this initializes once fo all, so if I use this same class twice it is not initialized twice with different inputs. I should rather create a function.
        #        class_experiment.__init__(self,filename,errfilename, modelfile,calplanck)
        self.filename=filename
        self.binningfile=nm.loadtxt((self.filename))[:,1]
        self.relcuts=nm.loadtxt((self.filename))[:,1]

        self.errfilename=errfilename
        self.model=nm.loadtxt(modelfile)
        #To check: now all the models for all data are rescaled by calplanck, and not just plik
        for ii in range(1,6):self.model[:,ii]=self.model[:,ii] #TBCHECKED
        self.lm={}
        self.lm['PP']=[8,2048]
        self.indmod={'TT':1,'TE':2,'EE':3,'PP':5}
        self.fsky={}

        lm= nm.loadtxt((self.filename))[:,3]
        fact=lm*(lm+1.)/nm.pi/2.
        temp=nm.loadtxt(self.filename)[:,4]
        err=nm.loadtxt(self.filename)[:,5]

        self.power=nm.transpose(nm.vstack((lm,temp,err)))
        print ('shape', self.power.shape)
        
        #        self.ind={'TT':1+0,'EE':1+1,'TE':1+3,'errpTT':1+0+6,'errpEE':1+1+6,'errpTE':1+3+6,'errnTT':1+0+6,'errnEE':1+1+6,'errnTE':1+3+6}


        self.ind={'PP':1,'errpPP':2,'errnPP':2}
        self.fsky['PP']=1
        self.binnedmodel=nm.loadtxt(binnedmodel)

        
# class class_lensing2(class_lensing):

#     def __init__(self,filename,errfilename, modelfile,calplanck):
#         class_lensing2.__init__(self,filename,errfilename, modelfile,calplanck)
        
