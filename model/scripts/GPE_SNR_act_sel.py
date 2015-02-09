#! Imports
import math
import numpy
import pylab
import os
import sys
import time

if len(sys.argv) != 1: mpiRun = True
else:                  mpiRun = False
start = time.time() 

 
# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])       
                   
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 
picture_dir=  '/'.join(os.getcwd().split('/')[0:-3]) + '/pictures'     
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
SPATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 

if mpiRun:
    LOAD=False
else:
    LOAD=False
SELECTIONTHR=2.  # Hz
NEURONMODELS=['SNR_aeif']
SYNAPSE_MODELS=['GPE_SNR_gaba_s_min', 
               'GPE_SNR_gaba_p']
#SYNAPSE_MODELS=[ 'GPE_SNR_gaba_s_min']

#def plot_example(GPE, SNR):
    
    
def plot_example_SNR(ax, SNR_list):
    time_bin=100
    
    colors = misc.make_N_colors('Blues', 5)
    colors=['g',colors[1]]   
    labels=['static','plastic']
    
    for color, label, SNR in zip(colors, labels, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax,
                          kwargs={'color':color})
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_xlim([2000,3500])
    #ax.set_ylim([10,30])
    #ax.text('Data 1')
    ax.text( 0.05, 0.45, 'static' , transform=ax.transAxes, **{ 'color' : colors[0] })  
    ax.text( 0.05, 0.25, 'plastic' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    #ax.text( 0.05, 0.25, 'Set 1' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    #ax.text( 0.05, 0.15, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 
    #ax.text( 0.05, 0.05, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[4] })     

def plot_example_raster_GPE(ax, GPE_list):
    time_bin=100
    ax_twinx=ax.my_twinx()
    GPE=GPE_list[0]

    GPE.signals['spikes'].raster_plot(id_list=GPE,
                                      display=ax_twinx,kwargs={'color':'k',
                                                               'zorder':1})  

    
    ax_twinx.set_ylim([0,10])

    ax_twinx.set_ylabel('Neuron id')
   # GPE.signals['spikes'].my_firing_rate( bin=time_bin, display=ax,
   #                                       kwargs={'color':'r',
   #                                               'linewidth':9,
   #                                               'zorder':20})
    GPE.signals['spikes'].my_firing_rate( bin=time_bin, display=ax,
                                          kwargs={'color':'k',
                                                  'linewidth':3,})
    ax.set_xlim([2000,3500])
    ax_twinx.set_xlim([2000,3500])
    ax.set_ylim([0,40])
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax_twinx.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency GPEs (Hz)')     
    ax.text( 0.7, 0.85, 'All' , transform=ax.transAxes, **{ 'color' : 'k' })  
    ax.text( 0.7, 0.75, 'Selected' , transform=ax.transAxes, **{ 'color' : 'grey' }) 
    
def plot_selection_vs_neurons1(ax, nbNeurons, meanRates):
    
    colors = misc.make_N_colors('Blues', 5)
    colors=['g','r', colors[1], colors[2], colors[3]]   
    
    ax.plot(nbNeurons,meanRates[0,:],**{'label':'Weak', 'color':colors[0]})  
    ax.plot(nbNeurons,meanRates[1,:],**{'label':'Strong', 'color':colors[1]})  
    ax.plot(nbNeurons,meanRates[2,:],**{'label':'Set 1', 'color':colors[2]})  
    ax.plot(nbNeurons,meanRates[3,:],**{'label':'Set 1+2', 'color':colors[3]})
    ax.plot(nbNeurons,meanRates[4,:],**{'label':'Set 2', 'color':colors[4]})
    ax.plot(nbNeurons,[SELECTIONTHR]*len(nbNeurons),'--k',**{'label':''})
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('GPEs (#)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_xlim([0,61])
    #ax.text('Data 1')
    ax.text( 0.7, 0.85, 'Weak' , transform=ax.transAxes, **{ 'color' : colors[0] })  
    ax.text( 0.7, 0.75, 'Strong' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.7, 0.65, 'Set 1' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.7, 0.55, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 
    ax.text( 0.7, 0.45, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[4] })  
    
def plot_selection_vs_neurons2(ax, nbNeurons, meanRates):
    
    colors = misc.make_N_colors('Blues', 5)
    colors=['g','r', colors[1], colors[2], colors[3]]   
    
    ax.plot(nbNeurons,meanRates[0,:],**{'label':'Weak', 'color':colors[0]})  
    ax.plot(nbNeurons,meanRates[1,:],**{'label':'Strong', 'color':colors[1]})  
    ax.plot(nbNeurons,meanRates[2,:],**{'label':'Set 1', 'color':colors[2]})  
    ax.plot(nbNeurons,meanRates[3,:],**{'label':'Set 1+2', 'color':colors[3]})
    ax.plot(nbNeurons,meanRates[4,:],**{'label':'Set 2', 'color':colors[4]})
    ax.plot(nbNeurons,[SELECTIONTHR]*len(nbNeurons),'--k',**{'label':''})
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('GPEs (#)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_xlim([0,61])

    ax.text( 0.7, 0.85, 'Weak' , transform=ax.transAxes, **{ 'color' : colors[0] })  
    ax.text( 0.7, 0.75, 'Strong' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.7, 0.65,'Set 1' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.7, 0.55, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 
    ax.text( 0.7, 0.45, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[4] })  

def plot_selection_vs_neurons_full(ax, hzs,data):
    colors = misc.make_N_colors('Blues', 5)
    colors=['g','r', colors[1], colors[2], colors[3]]   
    
    syn=SYNAPSE_MODELS
    ax.plot(hzs,data[syn[0]]['thrVec'],**{'label':'Weak', 'color':colors[0]})  
    ax.plot(hzs,data[syn[1]]['thrVec'],**{'label':'Strong', 'color':colors[1]})  
    ax.plot(hzs,data[syn[2]]['thrVec'],**{'label':'Set 1', 'color':colors[2]})  
    ax.plot(hzs,data[syn[3]]['thrVec'],**{'label':'Set 1+2', 'color':colors[3]})
    ax.plot(hzs,data[syn[4]]['thrVec'],**{'label':'Set 2', 'color':colors[4]})
    ax.set_xlabel('Frequency GPE (Hz)') 
    ax.set_ylabel('GPEs at thr (#)')
    ax.my_set_no_ticks( yticks=4, xticks = 5 ) 
    ax.set_xlim([0,51])
    ax.set_ylim([0,75])
    #ax.text('Data 1')
    ax.text( 0.7, 0.85, 'Weak' , transform=ax.transAxes, **{ 'color' : colors[0] })  
    ax.text( 0.7, 0.75, 'Strong' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.7, 0.65, 'Set 1' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.7, 0.55, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 
    ax.text( 0.7, 0.45, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[4] })  
    

def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURONMODELS )
    
    SNR = MyGroup( NEURONMODELS[0], 1, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    tb = ''     
    tb = tb + infoString
    
    tb = tb + '\n'

    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)

def simulate_example(hz=20):
    sname_nb=hz+1000
    
    saveAt=SPATH+'/'+SNAME+'-'+NEURONMODELS[0]+'-example.pkl'
    
    nGPE=100
    EXPERIMENTS=range(20)
    nSelected=50
    selectionRate=hz
    baseRate=20.
    model_list=models()
    I_e=-5.
    simTime=3500.
    selectionTime=500.
    selectionOnset=2500.
    
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( model_list, NEURONMODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS)       
 
    GPE_list=[] # GPE input for each experiment
    for iExp in EXPERIMENTS:
        GPE_list.append(MyPoissonInput( n=nGPE, sd=True, spath=SPATH, 
                                        sname_nb=sname_nb+iExp))
    
    SNR_list=[] # SNR groups for each synapse
    for iSyn, syn in enumerate(SYNAPSE_MODELS):
        SNR_list.append( MyGroup( NEURONMODELS[0], n=len(EXPERIMENTS), 
                       params={'I_e':I_e},
                       mm_dt = .1, record_from=[''],
                       spath=SPATH, sname_nb=sname_nb+iSyn) )
    if not LOAD:
        for iExp in EXPERIMENTS:
     
            GPE=GPE_list[iExp]
        
            # Base rate
            for id in GPE[:]:                 
                GPE.set_spike_times(id=id, rates=[baseRate], times=[1], 
                                t_stop=simTime)               
      
            # Selection        
            for id in GPE[:]: 
                GPE.set_spike_times(id=id, rates=[baseRate, selectionRate,
                                                  baseRate], 
                                    times=[1, selectionOnset,
                                           selectionTime+selectionOnset], 
                                    t_stop=simTime)     
        
                 
            for iSyn, syn in enumerate(SYNAPSE_MODELS):       
                    target=SNR_list[iSyn][iExp]
                    my_nest.ConvergentConnect(GPE[:], [target], model=syn)
                      
        my_nest.MySimulate( simTime )

    
        for GPE in GPE_list:
            GPE.get_signal( 's' ) # retrieve signal    
        
        for SNR in SNR_list:
            SNR.get_signal( 's' ) # retrieve signal 
        misc.pickle_save([GPE_list,SNR_list] , saveAt)

    if LOAD:
        GPE_list, SNR_list=misc.pickle_load(saveAt)
        
        
    s='\n'
    s=s+'Example:\n'
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( len(EXPERIMENTS) ),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'Base rate:', str ( baseRate ),  '#' )
    s = s + ' %s %5s %3s \n' % ( 'Selection rate:', str ( selectionRate ),  '#' )
      
    return GPE_list, SNR_list, s
    
def simulate_selection_vs_neurons(selRateInterval=[0.0,500.0], hz=20):    
    sname_nb=hz  
    
    nGPE=500
    nExp=5
    if hz>7:
        nMaxSelected=60
    else:
        nMaxSelected=100
    
    baseRate=0.1
    selectionRate=hz
    I_e=-5.
    
    simTime=3500.
    model_list=models()
    selectionTime=3000.
    selectionOnset=500.
    
       
    expParams=[]
    expIntervals=[]
    
    iSNR=0
    for syn in SYNAPSE_MODELS:
        for iSel in range(nMaxSelected):
            expIntervals.append([iSNR, iSNR+nExp])
            for iExp in range(nExp):
                expParams.append((syn, iSel, iExp, iSNR ))
                iSNR+=1
    
    synIntervals=[]
    iSNR=0
    for syn in SYNAPSE_MODELS:
        synIntervals.append([iSNR, iSNR+nMaxSelected])
        iSNR+=nMaxSelected  
    
    
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( model_list, NEURONMODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS)       
    

    SNR = MyGroup( NEURONMODELS[0], n=len(expParams), params={'I_e':I_e},
                        mm_dt = .1, record_from=[''],
                        spath=SPATH, sname_nb=sname_nb)
    
    sourceBack=[]
    sourceSel=[]
    for iExp in range(nExp):
        # Background 
        tmpSourceBack=[]
        for iGPE in range(nGPE-1):                  
            spikeTimes=misc.inh_poisson_spikes( [baseRate], [1],                        
                       t_stop=simTime, n_rep=nExp, seed=iGPE+10*iExp )    
                               
            if any(spikeTimes):  
                tmpSourceBack.extend( my_nest.Create('spike_generator', 
                                params={'spike_times':spikeTimes} ))
        sourceBack.append(tmpSourceBack)
 

        
    if not LOAD: 
        for syn, iSel, iExp, iSNR in expParams:       
            print 'Connect SNR '+ str(SNR[iSNR]) + ' ' + syn
            target=SNR[iSNR]
            my_nest.ConvergentConnect(sourceBack[iExp][0:nGPE-iSel], 
                                      [target], model=syn)
            my_nest.ConvergentConnect(sourceSel[iExp][0:iSel+1], 
                                      [target], model=syn)
                        
        my_nest.MySimulate( simTime )
            
        SNR.save_signal( 's' )
        SNR.get_signal( 's' ) # retrieve signal    
        
        #SNR.get_signal( 'v','V_m' ) # retrieve signal
        #SNR.signals['V_m'].plot()
        #SNR.signals['spikes'].raster_plot()
        #pylab.show()
                
    if LOAD:
        SNR.load_signal( 's' )
       
        
            #SNR.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
                
            #SNR.signals['V_m'].plot(id_list=[5])
            #SNR.['spikes'].raster_plot()
            #pylab.show()
    t1=selRateInterval[0]
    t2=selRateInterval[1]
    
    tmpMeanRates1=[]
    tmpMeanRates2=[]
    tmpMeanRates3=[]
    tmpMeanRates4=[]
    tmpMeanRates1=SNR.signals['spikes'].mean_rates(selectionOnset+t1, 
                                                   selectionOnset+t2)    
    for interval in expIntervals:
        tmpMeanRates3.append(numpy.mean(tmpMeanRates1[interval[0]:
                                                      interval[1]], axis=0))  
    
    for interval in synIntervals:          
        tmpMeanRates4.append(tmpMeanRates3[interval[0]:interval[1]])

    meanRates=numpy.array(tmpMeanRates4)
    nbNeurons=numpy.arange(1,nMaxSelected+1,1)
    
    s='\n'
    s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( nGPE ),  '#' )     
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( nExp ),  '#' )    
    s = s + ' %s %5s %3s \n' % ( 'Base rate:',   str ( baseRate),'Hz' )     
    s = s + ' %s %5s %3s \n' % ( 'Selection rate:', str ( selectionRate ), 'Hz' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( selectionTime ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'I_e:', str ( I_e ), 'pA' )
    

    infoString=s
    
    return nbNeurons, meanRates, infoString

def simulate_selection_vs_neurons_full(selRateInterval):
    
    s='-sel-'+ str(selRateInterval[0])+'-'+str(selRateInterval[1])
    saveAt=SPATH+'/'+SNAME+'-'+NEURONMODELS[0]+s+'.pkl'
    
    hzs=range(1,51)
    RETRIEVE=False
    
    if RETRIEVE:
        data={}
        for syn in SYNAPSE_MODELS:
            data[syn]={}
            data[syn]['rates']=[]
            data[syn]['selMat']=[]
            data[syn]['thrVec']=[]
            
            
        for hz in hzs:
            n, r, s = simulate_selection_vs_neurons(selRateInterval,hz)
            for r, syn in zip(r,SYNAPSE_MODELS):
                data[syn]['rates'].append(r)
        
        # Create matricies   
        for syn in data.keys():
            rates=data[syn]['rates']
            maxLen=0
            for r in rates:
                if len(r)>maxLen:
                    maxLen = len(r)
            for iR, r in  enumerate(rates):
                rates[iR]=numpy.append(r,numpy.zeros( (1,maxLen-len(r)) ))
                
         
        
        for syn in data.keys():
            r=numpy.array(data[syn]['rates'])
            selMat=numpy.array(r)
            thrVec=[]
            for i in range(r.shape[0]):
                p=True
                for j in range(r.shape[1]):
                    if SELECTIONTHR<r[i,j]:
                        selMat[i,j]=3
                    elif (SELECTIONTHR>=r[i,j]) and (SELECTIONTHR<r[i,j-1]) and p:
                        selMat[i,j]=2
                        thrVec.append(j+1)  # Neurons for threshold
                        p=False
                    else:
                        selMat[i,j]=1
                if p: 
                    thrVec.append(100)
    
                
            data[syn]['selMat']=selMat
            data[syn]['thrVec']=numpy.array(thrVec)
        
        if not mpiRun:
            misc.pickle_save(data, saveAt)
    
    elif not RETRIEVE:
        data=misc.pickle_load(saveAt)
    
    hzAt25actGPE=[]
    synNames1=''
    for syn in data.keys():
        data[syn]['thrVec']=numpy.array(data[syn]['thrVec'])
        n=data[syn]['thrVec'][data[syn]['thrVec']>25]
        hzAt25actGPE.append(len(n)+1)
        synNames1=synNames1+' '+syn[-3:]
    
    nbOfAt10hzGPE=[]
    synNames2=''
    for syn in data.keys():
        nbOfAt10hzGPE.append(data[syn]['thrVec'][9])
        synNames2=synNames2+' '+syn[-3:]
        
    s='\n'
    
    s = s + ' %s %5s %3s \n' % ( '25 act GPEs:\n%s\n'%(synNames1), str ( hzAt25actGPE ),  'Hz' )     
    s = s + ' %s %5s %3s \n' % ( '10 Hz  GPEs:\n%s\n'%(synNames2), str ( nbOfAt10hzGPE ),  '#' )     
    
    '''
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( nExp ),  '#' )    
    s = s + ' %s %5s %3s \n' % ( 'N max selected:', str ( nMaxSelected ),  '#' ) 
    s = s + ' %s %5s %3s \n' % ( 'Base rate:',   str ( baseRate),'Hz' )     
    s = s + ' %s %5s %3s \n' % ( 'Selection rate:', str ( selectionRate ), 'Hz' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( selectionTime ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( selectionOnset ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( simTime ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'I_e:', str ( I_e ), 'pA' )
    '''

    infoString=s
    
    return hzs, data, infoString
print 'Simulation'

# SIMULATION

# MOVIE  
'''      
for i, selRateInterval4 in enumerate([[200.0, 500.0],[500.0, 800.0],
                                      [800.0, 1100.0],[1100.0, 1400.0],
                                   [1400.0, 1700.0],[2000.0, 2300.0],
                                   [2300.0, 2600.0],[2600.0, 2900.0]]):
    hzs, data, s=simulate_selection_vs_neurons_full(selRateInterval4)
    
    fig = pylab.figure( facecolor = 'w' )
    ax=MyAxes(fig, [ .1, .1, .8,  .8 ] ) 
    plot_selection_vs_neurons_full(ax, hzs, data)
    ax.set_title('sel=%i-%i ms'%(selRateInterval4[0], selRateInterval4[1]), 
                 fontsize=12)
        
    filename = str('%03d' % i) + '.png'

    fig.savefig(filename, dpi=100)
#    fig.savefig( picture_dir + '/film/' +str('%03d' % i) + '.png', dpi = 100, format = 'png')

command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=4',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

os.spawnvp(os.P_WAIT, 'mencoder', command)
'''

infoString=''


GPE_list, SNR_list, s=simulate_example(hz=0)
infoString=infoString+s  

'''
selRateInterval4=[200.0, 500.0]
hzs, data, s=simulate_selection_vs_neurons_full(selRateInterval4)
infoString=infoString+s  

# simulate_selection_vs_neurons
LOAD=True
selRateInterval1=[0.0, 500.0]
nbNeurons1, meanRates1, s = simulate_selection_vs_neurons(selRateInterval1)

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )
'''
if not mpiRun:
    '''
    s = s + ' %s %5s %3s \n' % ( 'Sel thr:', str ( SELECTIONTHR ), 'Hz' )
    s = s + ' \n%s %5s %3s \n' % ( 'Sel interval 1:', str ( selRateInterval1 ), 'pA' )
        
    for mr, syn in zip(meanRates1, SYNAPSE_MODELS): 
       selNbNeruons=nbNeurons1[mr>SELECTIONTHR]
       s = s + ' %s %5s %3s \n' % ( syn+ 'sel' , 
                                    str ( selNbNeruons[-1]+1 ), '#' ) 
    infoString=infoString+s  
    
    
    selRateInterval2=[0.0, 200.0]
    nbNeurons2, meanRates2, s = simulate_selection_vs_neurons(selRateInterval2)
    s = ' \n%s %5s %3s \n' % ( 'Sel interval 2:', str ( selRateInterval2 ), 'pA' )       
    for mr, syn in zip(meanRates2, SYNAPSE_MODELS): 
       selNbNeruons=nbNeurons2[mr>SELECTIONTHR]
       s = s + ' %s %5s %3s \n' % ( syn+ 'sel' , 
                                    str ( selNbNeruons[-1]+1 ), '#' ) 
    infoString=infoString+s  
    
    
    
    selRateInterval3=[200.0, 500.0]
    nbNeurons3, meanRates3, s = simulate_selection_vs_neurons(selRateInterval3)
    s = ' \n%s %5s %3s \n' % ( 'Sel interval 3:', str ( selRateInterval3 ), 'pA' )       
    for mr, syn in zip(meanRates3, SYNAPSE_MODELS): 
       selNbNeruons=nbNeurons3[mr>SELECTIONTHR]
       s = s + ' %s %5s %3s \n' % ( syn+ 'sel' , 
                                    str ( selNbNeruons[-1]+1 ), '#' ) 
    infoString=infoString+s 
    '''
    
    # DISPLAY
    plot_settings.set_mode(mode='by_fontsize', w = 1100.0, h = 400.0, fontsize=12)
    font_size_text = 8
    fig = pylab.figure( facecolor = 'w' )
    
    ax_list = []
    ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
    ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .535, .6,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .535, .1,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 
              
    
    # Text
    ax=ax_list[0]
    plot_text(ax, infoString)
    
    ax=ax_list[1]
    plot_example_raster_GPE(ax, GPE_list)
    
    ax=ax_list[2]
    plot_example_SNR(ax, SNR_list)
    '''

    # plot_selection_vs_neurons
    ax=ax_list[3]
    plot_selection_vs_neurons1(ax, nbNeurons1, meanRates1)
    ax.set_title('sel=%i-%i ms'%(selRateInterval1[0], selRateInterval1[1]), 
                 fontsize=12)
    # plot_sele2tion_vs_neurons 
    ax=ax_list[4]   
    plot_selection_vs_neurons2(ax, nbNeurons2, meanRates2)  
    ax.set_title('sel=%i-%i ms'%(selRateInterval2[0], selRateInterval2[1]), 
                 fontsize=12)
    #ax.legend(numpoints=1, loc=[0.0,-1.5])  
    #plot_selection_vs_neurons_weak(ax, nbNeurons, meanRates)
    
    # plot_selection_vs_neurons 
    ax=ax_list[5]   
    plot_selection_vs_neurons1(ax, nbNeurons3, meanRates3)  
    ax.set_title('sel=%i-%i ms'%(selRateInterval3[0], selRateInterval3[1]), 
                 fontsize=12)
    
    ax=ax_list[6]  
    plot_selection_vs_neurons_full(ax, hzs, data)
    ax.set_title('sel=%i-%i ms'%(selRateInterval4[0], selRateInterval4[1]), 
                 fontsize=12)
    '''
    pylab.show()
    
    fig.savefig( picture_dir + '/' + SNAME  + '.svg', dpi = 500, format = 'svg')