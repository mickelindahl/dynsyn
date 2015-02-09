#!/usr/bin/python

# Module: plot_settings.py
# Author: Varun Hiremath <vh63@cornell.edu>
# Created: Thu,  2 Apr 2009 05:06:31 -0400

import pylab, math
from matplotlib.ticker import MaxNLocator
# Symbols
symbols = ['-','--','-.',':','.',',','o','^','v','<','>','s','+','x','D','d','1','2','3','4','h','H','p']
# Symbols + line
lps = [k+'-' for k in [',','.','o','^','v','<','>','s','+','x','D','d','1','2','3','4','h','H','p']]
# Colors
colors= ['b','g','r','c','m','y','k','w']




def push_remove_axis(ax, xaxis = 'push', yaxis = 'push', nb_xticks = 5, nb_yticks = 5):
    '''
    Push x and/or y axis  out from plot 10 points. Also remove upper and right
    frame.
    
    Inputs:
        ax        - axis handel
        xaxis     - 'push' or 'remove'
        yaxis     - 'push' or 'remove'
        nb_xticks - number of xticks to show
        nb_yticks - number of yticks to show 

    '''
    
  
    xlabel   = ax.get_xlabel()                                                  # Get xlabel
    ylabel   = ax.get_ylabel()                                                  # Get ylabel
    
    keep   = []
    push   = []
    remove = [ 'right', 'top' ]                                                 # Always remove up and right frame
    
    if   xaxis == 'keep': keep.append( 'bottom' )    
    elif xaxis == 'push': push.append( 'bottom' )
    else:                 remove.append( 'bottom' )
    
    if   yaxis == 'keep': keep.append( 'left' )   
    elif yaxis == 'push': push.append( 'left' )
    else:                 remove.append( 'left' )
    
    # Remove upper and right axes frames and move the bottom and left 
    for loc, spine in ax.spines.iteritems():
        if loc in push:     spine.set_position( ( 'outward', 10 ) )             # outward by 10 points
        elif loc in remove: spine.set_color( 'none' )                           # don't draw spine
        elif loc in keep: print 'Keeping axis'
        else:
            raise ValueError( 'unknown spine location: %s'%loc )
    
    if xaxis == 'remove':      
        ax.set_xticklabels( '' )                                                # turn of x ticks
        ax.xaxis.set_ticks_position( 'none' )                                   # turn of x tick labels
        ax.set_xlabel( '' )                                                     # turn of x label
    else: 
        ax.xaxis.set_ticks_position( 'bottom' )                                 # turn off x ticks where there is no spine
        ax.set_xlabel(xlabel,position=( 0.5, -0.2 ) )                           # add and reposition x label
    
    if yaxis == 'remove':
        ax.set_yticklabels( '' )                                                # turn of y ticks
        ax.yaxis.set_ticks_position( 'none' )                                   # turn of y tick labels
        ax.set_ylabel( '' )                                                     # turn of y label
    else: 
        ax.yaxis.set_ticks_position( 'left' )                                   # turn off y ticks where there is no spine
        ax.set_ylabel(ylabel,position=( -0.2, 0.5 ) )                           # Add and reposition y label
    

    ax.xaxis.set_major_locator( MaxNLocator( nb_xticks ) )                      # Set number of ticks x axis
    ax.yaxis.set_major_locator( MaxNLocator( nb_yticks ) )                      # Set number of ticks y axis
    
        

def get_figsize(fig_width_pt, w = None, h = None):
    if w == None and h == None:
        inches_per_pt = 1.0/72.0                # Convert pt to inch
        golden_mean = (math.sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*golden_mean      # height in inches
        fig_size =  [fig_width,fig_height]      # exact figsize
    else:
        inches_per_pt = 1.0/72.0                # Convert pt to inch
        fig_width  = w*inches_per_pt  # width in inches
        fig_height = h*inches_per_pt      # height in inches
        fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

# Publishable quality image settings for 2-column papers
params0 = {'backend': 'eps',
          'axes.labelsize': 6,
          'text.fontsize': 6,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'legend.pad': 0.1,    # empty space around the legend box
          'legend.fontsize': 5,
          'lines.markersize': 3,
          'font.size': 6,
          #'text.usetex': True,
          'figure.figsize': get_figsize(250),
          'figure.frameon':False}

# Medium sized images
params1 = {'backend': 'eps',
          'axes.labelsize': 8,
          'text.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'legend.pad': 0.1,     # empty space around the legend box
          'legend.fontsize': 8,
          'lines.markersize': 3,
          'font.size': 8,
          #'text.usetex': True,
          'figure.figsize': get_figsize(500),
          'figure.frameon':False}

# Large images (default)
params2 = {'backend': 'eps',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'legend.pad': 0.2,     # empty space around the legend box
          'legend.fontsize': 10,
           'lines.markersize': 3,
          'font.size': 10,
          #'text.usetex': True,
          'figure.figsize': get_figsize(800),
          'axes.frameon':False}

# My images (default)
params3 = {'backend': 'eps',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'legend.pad': 0.2,     # empty space around the legend box
          'legend.fontsize': 10,
          'lines.markersize': 3,
          'font.size': 10,
          #'text.usetex': True,
          'figure.figsize': get_figsize(800),
          'axes.frameon':False}

params4 = {'backend': 'eps',
          'axes.labelsize': 7,
          'text.fontsize': 7,
          'xtick.labelsize': 7,
          'ytick.labelsize': 7,
          'legend.pad': 0.1,     # empty space around the legend box
          'legend.fontsize': 7,
          'lines.markersize': 3,
          'font.size': 7,
          #'text.usetex': True,
          'figure.figsize': get_figsize(500),
          'figure.frameon':False}

params5 = {'backend': 'eps',
          'axes.labelsize': 15,
          'text.fontsize': 15,
          'title.fontsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15,
          'legend.pad': 0.2,    # empty space around the legend box
          'legend.fontsize': 15,
          'lines.markersize': 3,
          'font.size': 15,
          #'text.usetex': True,
          'figure.figsize': get_figsize(250),
          'figure.frameon':False}

params6 = {'backend': 'eps',
          'axes.labelsize': 20,
          'text.fontsize': 20,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'legend.pad': 0.2,     # empty space around the legend box
          'legend.fontsize': 20,
          'lines.markersize': 3,
          'font.size': 20,
          #'text.usetex': True,
          'figure.figsize': get_figsize(800),
          'axes.frameon':False}

def set_mode(mode='large', pl=pylab, w = None, h = None):
    if mode == 'my_large':
        params6['figure.figsize'] = get_figsize( 800, w = w, h = h)
        pl.rcParams.update(params6)
    if mode == 'my_publish':
        params3['figure.figsize'] = get_figsize( 800, w = w, h = h)
        pl.rcParams.update(params3)
    elif mode == 'my_medium':
        params5['figure.figsize'] = get_figsize( 500, w = w, h = h)
        pl.rcParams.update(params5)
    elif mode == 'my_low':
        params4['figure.figsize'] = get_figsize( 250, w = w, h = h)
        pl.rcParams.update(params4)    
    elif mode == "publish":
        pl.rcParams.update(params0)
    elif mode == "medium":
        pl.rcParams.update(params1)
    #else:
    # pl.rcParams.update(params2)

def set_figsize(fig_width_pt):
    pylab.rcParams['figure.figsize'] = get_figsize(fig_width_pt)





