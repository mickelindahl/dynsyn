import numpy
import pylab
import os
import sys
import time as ttime

# Get directory where model and code resides 
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])    
code_dir=  '/'.join(os.getcwd().split('/')[0:-2])  

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 

from simulation_utils import simulate_MSN
from src import misc

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
RUN=True

# 48 minuter 200x10x10

HZS=numpy.linspace(5,50,10)
N_MSNS=numpy.linspace(5,50,10)
t=ttime.time()
save_result_at=OUTPUT_PATH+'/main_mutual_information_raw_data.pkl'
if RUN is False:
    count_dic={}
    for i_syn, syn in enumerate(['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_mid',  'MSN_SNR_gaba_s_max', 
                 'MSN_SNR_gaba_p1']):
        count_dic[i_syn]={}
        
        for hz in numpy.linspace(5,50,10):
            count_dic[i_syn][hz]={}
            for N_MSN in numpy.linspace(5,50,10):
                count_dic[i_syn][hz][N_MSN]=[]
                for i in range(200):
                    c, time= simulate_MSN(int(N_MSN), ['SNR_izh'], [syn], 
                                          sim_time=1000, burst_onset=700, burst_offset=1000, 
                                          burst_rate=hz, threads=1)
                    count_dic[i_syn][hz][N_MSN].extend(c)
                count_dic[i_syn][hz][N_MSN]=numpy.array(count_dic[i_syn][hz][N_MSN])
    
    
    
    misc.pickle_save([count_dic, time],save_result_at)

else:
    #count_dic, time= misc.pickle_load(save_result_at)
    pass

save_result_at=OUTPUT_PATH+'/main_mutual_information_prob.pkl'
if RUN is False:
    c_prob=[]
    c_sum=[]
    for i in sorted(list(count_dic.keys())):
        c_prob.append([])
        c_sum.append([])
        for j, hz in enumerate(sorted(list(count_dic[i].keys()))):
            c_prob[i].append([])
            c_sum[i].append([])
            for N in sorted(list(count_dic[i][hz].keys())):
                
                c=count_dic[i][hz][N][2::3,:] # Pick out results for SNr
                c_conv=misc.convolve(c,bin_extent=100, kernel_type='rectangle')
                c_prob[i][j].append(numpy.sum(c_conv<1,axis=0)/float(c_conv.shape[0]))
                c_sum[i][j].append(numpy.sum(c,axis=0))
                
    c_prob=numpy.array(c_prob)
    c_sum=numpy.array(c_sum)

    misc.pickle_save([c_prob, c_sum, time],save_result_at)
else:
    c_prob, c_sum, time= misc.pickle_load(save_result_at)            




from numpy import *
import pylab as p
#import matplotlib.axes3d as p3
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm

#ax2=pylab.subplot(111)

for i in range(c_sum.shape[0]):

    #
    x=[]
    y=[]
    z=[]
    #fig = pylab.figure()
    #ax1=pylab.subplot(111)
    for j in range(c_sum.shape[1]):
        x.append([])
        y.append([])
        z.append([])
        for k in range(c_sum.shape[2]):
            x[j].append(HZS[j])
            y[j].append(N_MSNS[k])
            z[j].append(c_sum[i,j,k,770])
            #z[j].append(c_prob[i,j,k,770])
           
            if (i==1) and (j==1) and (k==9):  
                fig = pylab.figure() 
                ax=pylab.subplot(111)
                ax.plot(time, c_sum[i,j,k,:])
                ax.plot(time, misc.convolve(c_sum[i,j,k,:], bin_extent=100, kernel_type='rectangle', single=True).transpose())
                pylab.show()
                ax.set_xlabel('Spike count')
            #ax1.plot(c_sum[i,j,k,:])
            #ax1.plot(c_prob[i,j,k,:])
    x=numpy.array(x)
    y=numpy.array(y)
    z=numpy.array(z)   
#    fig = pylab.figure()     
#    ax = p3.Axes3D(fig)
#    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet)
#    ax.view_init(50, 135)
#    
    
    fig = pylab.figure() 
    ax=pylab.subplot(111)
    for k in range(c_sum.shape[1]):
        ax.plot(x[:,k],z[:,k])
    ax.set_xlabel('Firing rate (spikes/s)')
    
    fig = pylab.figure() 
    ax=pylab.subplot(111)
    for k in range(c_sum.shape[2]):
        ax.plot(y[k,:],z[k,:])
    ax.set_xlabel('Number of MSNs')
    
    #ax1.pcolor(numpy.array(z))
    #pylab.colorbar()

pylab.show()

count=numpy.array(count_dic[1][5][5])
save_result_at=OUTPUT_PATH+'/main_mutual_information_convolved_raw_data.pkl'   
if RUN is True:
    count_conv=[]
    count_conv=misc.convolve(count, bin_extent=10, kernel_type='rectangle')
  
    misc.pickle_save([count_conv],save_result_at)
else:
    count_conv= misc.pickle_load(save_result_at)
count=count_conv
    
save_result_at=OUTPUT_PATH+'/main_mutual_information_count_sr.pkl'   
if RUN is True:        
    data=[]
    min_data=[]
    max_data=[]
        
    for i in range(3):
        data.append(numpy.array(count[i::3,:]))
        min_data.append(numpy.floor(numpy.min(count[i::3,:])))
        max_data.append(numpy.ceil(numpy.max(count[i::3,:])))
    #data[1]=data[0]+data[1]   
    
    count_rs_list=[]    

    print max_data
    for k in range(len(time)):
        count_rs=numpy.zeros((max_data[1]+1,max_data[2]+1))
        
        nx=range(min_data[1], max_data[1]+2)
        ny=range(min_data[2], max_data[2]+2)
        
        count_rs, xedges, yedges = numpy.histogram2d(data[1][:,k], data[2][:,k], bins=(nx, ny))
     
            
        count_rs_list.append(count_rs)  
    misc.pickle_save([count_rs_list,data], save_result_at)
        
else:
    count_rs_list, data=misc.pickle_load(save_result_at)


ff=misc.fano_factor(data[2])
a=numpy.ones((2,2))
b=numpy.ones((3,3))
c=numpy.array([[1,2,3],[2,1,3],[3,2,1]])
d=numpy.array([[2,1],[2,2]])
mi_test=misc.mutual_information([d,c, a,b])  
    
mi=misc.mutual_information(count_rs_list)            
    
max_data=[]
for i in range(3):
    data.append(numpy.array(count[i::3]))
    max_data.append(numpy.max(count[i::3]))
    
print count_rs_list[0].shape, max_data, len(count_rs_list)
    
pylab.figure()
pylab.plot(time, ff) 


pylab.figure()
for i in range(3):
    pylab.subplot(3,1,i+1)
    pylab.plot(time, numpy.sum(data[i], axis=0)) 
    
print 'Simulation sime', (t-ttime.time()) / 60., 'min'
pylab.show()



