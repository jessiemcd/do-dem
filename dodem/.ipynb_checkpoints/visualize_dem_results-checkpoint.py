import pickle
from matplotlib import pyplot as plt
import numpy as np
import copy

import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter, ScalarFormatter
from astropy import units as u
import os

from scipy import optimize


def make_timestring(time):
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')

    
    return timestring+'_'+stopstring


def load_DEM(filename):
    """
    Take in a DEM file (saved output from dodem.dodem) + return dictionary of results + inputs.
    
    CURRENT FUNCTIONALITY: Enter name of file (pickle) with DEM results.
    
    OLD FUNCTIONALITY: Instead of file name, requests time interval and fpm (standardized file 
                        names+locations assumed).
    """
    
    try:
        with open(filename, 'rb') as f:
            data = pickle.load(f)
    except FileNotFoundError:
        print('No DEM file: ', filename)
        return

    time = data['time_interval']
    timestring=make_timestring(time)
    
    return data, timestring, time

def plot_DEM(data, title='', plotMK=False, fill_color='lightcoral'):
    """
    Takes in DEM-output dictionary, and makes a nice plot.
    Note: to plot two DEMs at once for comparison, see compare_DEMs() below. 
    
    Doesn't return anything, but saves the plot (& you can enter a title). 
    
    """
    
    #Turn off for linear y scale (not a very pretty option)
    logplot=True

    fig = plt.figure(figsize=(18, 10), tight_layout = {'pad': 1})

    fig.add_subplot(1,2,1)

    clrs = comparison_colors(data)
    edem=np.array(data['edem'])
    dem=np.array(data['DEM'])
    
    #Plot loci curves for each instrument:
    for i in np.arange(len(data['chanax'])):
        if plotMK == True:
            plt.loglog(data['ts_'],data['dn_in'][i]/data['trmatrix'][:,i],
                       label=data['chanax'][i],color=clrs[i],lw=4)
        else:
            if logplot:
                plt.semilogy(data['ts_'],data['dn_in'][i]/data['trmatrix'][:,i],
                         label=data['chanax'][i],color=clrs[i],lw=4)
            else:
                plt.plot(data['ts_'],data1['dn_in'][i]/data['trmatrix'][:,i],
                         label=data['chanax'][i],color=clrs[i],lw=4)
                

    if data['chisq'].ndim == 1:
        chi_sq = data['chisq'][0]
    else:
        chi_sq = data['chisq']
                
    #Plot DEMs (result: plot markers. uncertainty range: shaded region)        
    if len(edem.shape) == 2:
        #Asymmetric error
        plt.fill_between(data['ts'], dem-edem[0,:], dem+edem[1,:], color=fill_color, alpha=0.5)
    if len(edem.shape) == 1:
        #Symmetric error 
        plt.fill_between(data['ts'], dem-edem, dem+edem, color=fill_color, alpha=0.5)
        
    if 'xdem_error' in data:
        plt.errorbar(data['ts'],dem,xerr=data['xdem_error'],fmt='^r',label='AXN $\chi^2 =$ {:0.2f}'.format(chi_sq))
    else:    
        plt.scatter(data['ts'],dem,marker='^',label='AXN $\chi^2 =$ {:0.2f}'.format(chi_sq), color='red')
    
    plt.grid(True,which='both',lw=0.5,color='gainsboro')
    if logplot:
        plt.ylim([1e20,9e29])
    else:
        plt.ylim([1e17,1e27])
    plt.xlim([5.5,8])
    if plotMK == True:
        plt.xlim([min(data['ts']),9])
    else:
        plt.xlim([min(data['ts']),max(data['ts'])])

    plt.title('DEM With Instrument Loci Curves', fontsize = 30)
    plt.legend(ncol=3, fontsize = 15)
    if plotMK == True:
        plt.xlabel('T [MK]', fontsize=20)
    else:
        plt.xlabel('$\mathrm{\log_{10}T\;[K]}$', fontsize=20)
    plt.ylabel('$Emission\;Measure\;[cm^{-5}]$', fontsize=20)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)

    ax = fig.add_subplot(1,2,2)
    
    residuals = data['dn_reg']/data['dn_in']
    eresiduals = data['edn']
    nf = len(residuals)

    plt.errorbar(np.arange(0,nf,1),residuals, yerr=eresiduals, fmt='^r',
                     label='$\chi^2 =$ {:0.2f}'.format(chi_sq)+', '+data['edn_string'])
    
    #Line at y=1 (ideal)
    plt.plot([-1,nf+1],[1,1],'--',color='grey')


    ax.set_ylim([0.5,2])
    ax.set_xlim([-0.5,nf])
    ax.set_yscale('log')
    ax.set_xticks(np.arange(nf),data['chanax'],rotation=90, fontsize=15)
    ax.set_yticks([0.5,1,2], ['0.5','1','2'], fontsize=15)
    ax.set_xlabel('Channel', fontsize=20)
    ax.set_ylabel('DN$_\mathrm{(DEM\ predicted)}$/DN$_\mathrm{(measured)}$', fontsize=20)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize = 15)
    
    plt.savefig(title+'_DEM_plot.png')
    plt.close(fig)
    
    

def compare_DEMs(data1, data2, timestring1, timestring2, title1='', title2='', plotMK=False, plot=True):
    """
    Applicability:
    ---------------
            WORKS if the two DEMs are over the SAME temperature array (i.e. to compare two same-method
            DEMs at different times).

            ALSO WORKS if they have the same lower temperature bound and the same difference between the 
            first two temperatures (assuming this means they have the same uniform temperature step), even
            if one array contains more higher temperatures. In this case, we pad the shorter DEM result with 
            zeros for plotting, etc. (This is to compare the DEMs over the same time range, but different 
            temperature ranges.)

            FAILS in any other situation.
    
    Allows you to "name" the two DEMs via title keywords.
    
    """
    tempcomp=False

    if np.array_equal(data1['ts'], data2['ts']) == False:
        tempcomp=True
        print('DEMs are done over different temperature arrays.')
        if data1['ts'][0] == data2['ts'][0] and data1['ts'][1] == data2['ts'][1]:
            print('Same low-temp bound and step size.')
            diff = data1['ts'][1]-data1['ts'][0]
            lendiff = len(data1['ts'])-len(data2['ts'])
            if lendiff > 0:
                print('data1 is over a larger temperature range.')
                data2['ts'] = np.append(data2['ts'], data1['ts'][-lendiff:])
                if np.array_equal(data1['ts'],data2['ts']):
                    data2['DEM'] = np.append(data2['DEM'], np.zeros(lendiff))
                    data2['edem'][0] = np.append(data2['edem'][0], np.zeros(lendiff))
                    data2['edem'][1] = np.append(data2['edem'][1], np.zeros(lendiff))
                else:
                    print('data2 is over a larger temperature range.')
                    lendiff=abs(lendiff)
                    data1['ts'] = np.append(data1['ts'], data2['ts'][-lendiff:])
                    if np.array_equal(data1['ts'],data2['ts']):
                        data1['DEM'] = np.append(data1['DEM'], np.zeros(lendiff))
                        data1['edem'][0] = np.append(data1['edem'][0], np.zeros(lendiff))
                        data1['edem'][1] = np.append(data1['edem'][1], np.zeros(lendiff))
            
        else:
            print("Results DO NOT have the same low-temp bound and step size.")
            print("Quitting")
            return
    
    clrs1 = comparison_colors(data1)
    clrs2 = comparison_colors(data2)
    edem1=np.array(data1['edem'])
    edem2=np.array(data2['edem'])
    dem1=np.array(data1['DEM'])
    dem2=np.array(data2['DEM'])
    
    
    #Locations where the two DEMs are consistent within mutual uncertainties
    consistent = np.zeros(len(dem1))
    for i in range(0, len(dem1)):
        evals1 = [dem1[i]-edem1[:,i][0], dem1[i]+edem1[:,i][1]]
        evals2 = [dem2[i]-edem2[:,i][0], dem2[i]+edem2[:,i][1]]
        
        if evals2[0] < evals1[0] < evals2[1] or evals1[0] < evals2[0] < evals1[1]:
            consistent[i] = 1
            
    #Make DEM Comparison plot
    if plot:
        
        logplot=True
        
        fig = plt.figure(figsize=(25, 12), tight_layout = {'pad': 1})

        fig.add_subplot(1,2,1)

        #Vertical lines where DEMs are consistent
        #for c in range(0, len(consistent)):
        #    if consistent[c]==1:
        #        plt.axvline(data1['ts'][c], color='Black', linestyle='dashed')
                #print(data1['ts'][c])
        
        #Plot loci curves for each instrument, DEM 1:
        for i in np.arange(len(data1['chanax'])):
            if plotMK == True:
                plt.loglog(data1['ts_'],data1['dn_in'][i]/data1['trmatrix'][:,i],
                           label=data1['chanax'][i],color=clrs1[i],lw=4)
                            #+'-'+str(1)
            else:
                if logplot:
                    plt.semilogy(data1['ts_'],data1['dn_in'][i]/data1['trmatrix'][:,i],
                             label=data1['chanax'][i],color=clrs1[i],lw=4)
                else:
                    plt.plot(data1['ts_'],data1['dn_in'][i]/data1['trmatrix'][:,i],
                             label=data1['chanax'][i],color=clrs1[i],lw=4)
                    
        #Plot loci curves for each instrument, DEM 2:
        for i in np.arange(len(data2['chanax'])):
            if plotMK == True:
                plt.loglog(data2['ts_'],data2['dn_in'][i]/data2['trmatrix'][:,i],
                           #label=data2['chanax'][i]+'-'+str(2), 
                           color=clrs2[i],lw=4, linestyle='dashed')
            else:
                if logplot:
                    plt.semilogy(data2['ts_'],data2['dn_in'][i]/data2['trmatrix'][:,i],
                             #label=data2['chanax'][i]+'-'+str(2), 
                                 color=clrs2[i],lw=4, linestyle='dashed')
                else:
                    plt.plot(data2['ts_'],data2['dn_in'][i]/data2['trmatrix'][:,i],
                             #label=data2['chanax'][i]+'-'+str(2), 
                             color=clrs2[i],lw=4, linestyle='dashed')

                    
        if data1['chisq'].ndim == 1:
            chi_sq1 = data1['chisq'][0]
        else:
            chi_sq1 = data1['chisq']

        if data2['chisq'].ndim == 1:
            chi_sq2 = data2['chisq'][0]
        else:
            chi_sq2 = data2['chisq']
        
        #Check if the DEMs have their own defined colors (each method has a color to distinguish plots)
        fill1='lightcoral'
        fill2='lightblue'
        if 'fill_color' in data1:
            fill1 = data1['fill_color']
        if 'fill_color' in data2:
            fill2 = data2['fill_color']
          
        #print(fill1, fill2)
        #If both of the DEMs had the same method, change one color to tell them apart still!

        if fill1==fill2:
            #print('hi')
            fill2='lightblue'

            
        if 'edn_string' in data1:
            edn_string1 = data1['edn_string']
        if 'edn_string' in data2:
            edn_string2 = data2['edn_string']
            
        #Plot DEMs (result: plot markers. uncertainty range: shaded region)
        plt.fill_between(data1['ts'], dem1-edem1[0,:], dem1+edem1[1,:], color=fill1, alpha=0.5)
        plt.scatter(data1['ts'],dem1,marker='^',
                    label='1: $\chi^2 =$ {:0.2f}'.format(chi_sq1), color=fill1)
                    #+'-'+str(1)

        #color=fill2,alpha=0.5,    
        plt.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2, hatch='/')
        plt.scatter(data1['ts'],dem2,marker='^',
                    label='2: $\chi^2 =$ {:0.2f}'.format(chi_sq2), color=fill2)

        plt.grid(True,which='both',lw=0.5,color='gainsboro')
        if logplot:
            plt.ylim([1e20,9e29])
        else:
            plt.ylim([1e17,1e27])
        plt.xlim([5.5,8])
        if plotMK == True:
            plt.xlim([min(data1['ts']),9])
        else:
            plt.xlim([min(data1['ts']),max(data1['ts'])])

        plt.title('Compare DEMs - 1: '+title1+', 2: '+title2, fontsize = 30)
        plt.legend(ncol=2, fontsize = 20, loc='lower left')
        if plotMK == True:
            plt.xlabel('T [MK]', fontsize=20)
        else:
            plt.xlabel('$\mathrm{\log_{10}T\;[K]}$', fontsize=30)
        plt.ylabel('$Emission\;Measure\;[cm^{-5}]$', fontsize=30)
        plt.yticks(fontsize=25)
        plt.xticks(fontsize=25)

        ax = fig.add_subplot(1,2,2)
        
        #Get residuals, ready for comparison plot.
        residuals1, indices1, residuals2, indices2, inst = comparison_instruments(data1, data2)
        
        nf = len(inst)
        
        indices1 = [i-0.05 for i in indices1]
        indices2 = [i+0.05 for i in indices2]
        
        ax.errorbar(indices1,residuals1, yerr=data1['edn'], fmt='^', color=fill1,
                     label='1: $\chi^2 =$ {:0.2f}'.format(chi_sq1) )#+' '+edn_string1)
        
        ax.errorbar(indices2,residuals2, yerr=data2['edn'], fmt='^', color=fill2,
                     label='2: $\chi^2 =$ {:0.2f}'.format(chi_sq2) )#+' '+edn_string2)
        
        #Line at y=1 (ideal)
        ax.plot([-1,nf+1],[1,1],'--',color='grey')

        ax.set_ylim([0.5,2])
        ax.set_xlim([-0.5,nf])
        ax.set_yscale('log')
        ax.set_title('DEM Residuals', fontsize=30)
        ax.set_xticks(np.arange(nf),inst,rotation=45, fontsize=25, ha='right')
        ax.tick_params(axis='x', length=20)
        ax.set_yticks([0.5,1,2], ['0.5','1','2'], fontsize=25)
        #ax.set_xlabel('Channel', fontsize=30)
        ax.set_ylabel('DN$_\mathrm{(DEM\ predicted)}$/DN$_\mathrm{(measured)}$', fontsize=30)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.legend(fontsize = 20)
        
        if tempcomp:
            plt.savefig('./'+timestring1+'/DEM_comparison_'+title1+'_vs_'+title2+'temp_range_comp.png')
            plt.close(fig)
        else:
            plt.savefig('DEM_comparison_'+title1+'_vs_'+title2+'.png')
            
#     m1 = DEMmax(data1['ts'], data1['DEM'], wind=3, plot=False)
#     m2 = DEMmax(data2['ts'], data2['DEM'], wind=3, plot=False)

#     print('DEM 1 is a maximum at: log(T)=', m1, 'OR, ', 10**m1/1e6, ' MK')
#     print('DEM 2 is a maximum at: log(T)=', m2, 'OR, ', 10**m2/1e6, ' MK')
        
    return consistent




def compare_DEMs_4panel(data1, data2, timestring1, timestring2, title1='', title2='', plotMK=False, plot=True,
                       fill1='lightcoral', fill2='lightskyblue', fixfill=False, peak_inset=False, high_inset=False):
    """
    Applicability:
    ---------------
            WORKS if the two DEMs are over the SAME temperature array (i.e. to compare two same-method
            DEMs at different times).

            ALSO WORKS if they have the same lower temperature bound and the same difference between the 
            first two temperatures (assuming this means they have the same uniform temperature step), even
            if one array contains more higher temperatures. In this case, we pad the shorter DEM result with 
            zeros for plotting, etc. (This is to compare the DEMs over the same time range, but different 
            temperature ranges.)

            FAILS in any other situation.
    
    Allows you to "name" the two DEMs via title keywords.
    
    """
    tempcomp=False

    if np.array_equal(data1['ts'], data2['ts']) == False:
        tempcomp=True
        print('DEMs are done over different temperature arrays.')
        if data1['ts'][0] == data2['ts'][0] and data1['ts'][1] == data2['ts'][1]:
            print('Same low-temp bound and step size.')
            diff = data1['ts'][1]-data1['ts'][0]
            lendiff = len(data1['ts'])-len(data2['ts'])
            if lendiff > 0:
                print('data1 is over a larger temperature range.')
                data2['ts'] = np.append(data2['ts'], data1['ts'][-lendiff:])
                if np.array_equal(data1['ts'],data2['ts']):
                    data2['DEM'] = np.append(data2['DEM'], np.zeros(lendiff))
                    data2['edem'][0] = np.append(data2['edem'][0], np.zeros(lendiff))
                    data2['edem'][1] = np.append(data2['edem'][1], np.zeros(lendiff))
                else:
                    print('data2 is over a larger temperature range.')
                    lendiff=abs(lendiff)
                    data1['ts'] = np.append(data1['ts'], data2['ts'][-lendiff:])
                    if np.array_equal(data1['ts'],data2['ts']):
                        data1['DEM'] = np.append(data1['DEM'], np.zeros(lendiff))
                        data1['edem'][0] = np.append(data1['edem'][0], np.zeros(lendiff))
                        data1['edem'][1] = np.append(data1['edem'][1], np.zeros(lendiff))
            
        else:
            print("Results DO NOT have the same low-temp bound and step size.")
            print("Quitting")
            return
    
    clrs1 = comparison_colors(data1)
    clrs2 = comparison_colors(data2)
    edem1=np.array(data1['edem'])
    edem2=np.array(data2['edem'])
    dem1=np.array(data1['DEM'])
    dem2=np.array(data2['DEM'])
    
    
    #Locations where the two DEMs are consistent within mutual uncertainties
    consistent = np.zeros(len(dem1))
    for i in range(0, len(dem1)):
        evals1 = [dem1[i]-edem1[:,i][0], dem1[i]+edem1[:,i][1]]
        evals2 = [dem2[i]-edem2[:,i][0], dem2[i]+edem2[:,i][1]]
        
        if evals2[0] < evals1[0] < evals2[1] or evals1[0] < evals2[0] < evals1[1]:
            consistent[i] = 1
            
    #Make DEM Comparison plot
    if plot:
        
        logplot=True
        
        fig, axes = plt.subplots(2,2, figsize=(25, 18), tight_layout = {'pad': 1})

        ax=axes[1,0]
        ax1=axes[1,1]
        ax2=axes[0,0]
        ax3=axes[0,1]

        #Vertical lines where DEMs are consistent
        #for c in range(0, len(consistent)):
        #    if consistent[c]==1:
        #        plt.axvline(data1['ts'][c], color='Black', linestyle='dashed')
                #print(data1['ts'][c])
        
        #Plot loci curves for each instrument, DEM 1:
        for i in np.arange(len(data1['chanax'])):
            if plotMK == True:
                ax2.loglog(data1['ts_'],data1['dn_in'][i]/data1['trmatrix'][:,i],
                           label=data1['chanax'][i],color=clrs1[i],lw=4)
                            #+'-'+str(1)
            else:
                if logplot:
                    ax2.semilogy(data1['ts_'],data1['dn_in'][i]/data1['trmatrix'][:,i],
                             label=data1['chanax'][i],color=clrs1[i],lw=4)
                else:
                    ax2.plot(data1['ts_'],data1['dn_in'][i]/data1['trmatrix'][:,i],
                             label=data1['chanax'][i],color=clrs1[i],lw=4)
                    
        #Plot loci curves for each instrument, DEM 2:
        for i in np.arange(len(data2['chanax'])):
            if plotMK == True:
                ax3.loglog(data2['ts_'],data2['dn_in'][i]/data2['trmatrix'][:,i],
                           label=data2['chanax'][i], #+'-'+str(2), 
                           color=clrs2[i],lw=4)
            else:
                if logplot:
                    ax3.semilogy(data2['ts_'],data2['dn_in'][i]/data2['trmatrix'][:,i],
                             label=data2['chanax'][i], #+'-'+str(2), 
                                 color=clrs2[i],lw=4)
                else:
                    ax3.plot(data2['ts_'],data2['dn_in'][i]/data2['trmatrix'][:,i],
                             label=data2['chanax'][i], #+'-'+str(2), 
                             color=clrs2[i],lw=4)

                    
        if data1['chisq'].ndim == 1:
            chi_sq1 = data1['chisq'][0]
        else:
            chi_sq1 = data1['chisq']

        if data2['chisq'].ndim == 1:
            chi_sq2 = data2['chisq'][0]
        else:
            chi_sq2 = data2['chisq']
        
        if fixfill==False:
            #Check if the DEMs have their own defined colors (each method has a color to distinguish plots)
            if 'fill_color' in data1:
                fill1 = data1['fill_color']
            if 'fill_color' in data2:
                fill2 = data2['fill_color']

            #print(fill1, fill2)
            #If both of the DEMs had the same method, change one color to tell them apart still!

            if fill1==fill2:
                #print('hi')
                fill2='lightskyblue'

            
        if 'edn_string' in data1:
            edn_string1 = data1['edn_string']
        if 'edn_string' in data2:
            edn_string2 = data2['edn_string']
            
        #Plot DEMs (result: plot markers. uncertainty range: shaded region)
        ax.fill_between(data1['ts'], dem1-edem1[0,:], dem1+edem1[1,:], color=fill1, alpha=0.5)
        ax.scatter(data1['ts'],dem1,marker='^',
                    label='1: '+title1+' ($\chi^2 =$ {:0.2f})'.format(chi_sq1), color=fill1, edgecolors= "black")
        ax2.fill_between(data1['ts'], dem1-edem1[0,:], dem1+edem1[1,:], color=fill1, alpha=0.5)
        ax2.scatter(data1['ts'],dem1,marker='^',
                    label='1: $\chi^2 =$ {:0.2f}'.format(chi_sq1), color=fill1)

        #color=fill2,alpha=0.5,    
        ax.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2)#, hatch='/')
        ax.scatter(data1['ts'],dem2,marker='^',
                    label='2: '+title2+' ($\chi^2 =$ {:0.2f})'.format(chi_sq2), color=fill2, edgecolors= "black")
        ax3.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2)#, hatch='/')
        ax3.scatter(data1['ts'],dem2,marker='^',
                    label='2: $\chi^2 =$ {:0.2f}'.format(chi_sq2), color=fill2)
        
        
        if peak_inset:
            x1, x2, y1, y2 = 6.35, 6.65, 1e25, 1e27  # subregion of the original image - little box around peak
            #x1, x2, y1, y2 = 6.35, 6.8, 1e24, 1e27
            axins = ax.inset_axes(
                [5.7, 2e21, 0.6, 2e25], #Little box around peak
                #[5.7, 2e21, 0.675, 5e25], #Little box around peak
                xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[], transform=ax.transData)
            axins.fill_between(data1['ts'], dem1-edem1[0,:], dem1+edem1[1,:], color=fill1, alpha=0.5)
            axins.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2)#, hatch='/')
            axins.scatter(data1['ts'],dem1,marker='^',
                    label='1: '+title1+' ($\chi^2 =$ {:0.2f})'.format(chi_sq1), color=fill1, edgecolors= "black")
            #axins.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2)#, hatch='/')
            axins.scatter(data1['ts'],dem2,marker='^',
                    label='2: '+title2+' ($\chi^2 =$ {:0.2f})'.format(chi_sq2), color=fill2, edgecolors= "black")
            axins.set_yscale('log')
            axins.yaxis.set_major_formatter(NullFormatter())
            ax.indicate_inset_zoom(axins, edgecolor="black")

        if high_inset:
            x1, x2, y1, y2 = 6.8, 7.2, 1e22, 1e25  # subregion of the original image - little box around peak
            #x1, x2, y1, y2 = 6.35, 6.8, 1e24, 1e27
            axins = ax.inset_axes(
                [6.0, 2e21, 0.6, 5e25], #Little box around peak
                #[5.7, 2e21, 0.675, 5e25], #Little box around peak
                xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[], transform=ax.transData)
            axins.fill_between(data1['ts'], dem1-edem1[0,:], dem1+edem1[1,:], color=fill1, alpha=0.5)
            axins.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2)#, hatch='/')
            axins.scatter(data1['ts'],dem1,marker='^',
                    label='1: '+title1+' ($\chi^2 =$ {:0.2f})'.format(chi_sq1), color=fill1, edgecolors= "black")
            #axins.fill_between(data1['ts'], dem2-edem2[0,:], dem2+edem2[1,:], alpha=0.5, color=fill2)#, hatch='/')
            axins.scatter(data1['ts'],dem2,marker='^',
                    label='2: '+title2+' ($\chi^2 =$ {:0.2f})'.format(chi_sq2), color=fill2, edgecolors= "black")
            axins.set_yscale('log')
            axins.yaxis.set_major_formatter(NullFormatter())
            ax.indicate_inset_zoom(axins, edgecolor="black")            
        
        
        ax.set_yscale('log')
        ax.grid(True,which='both',lw=0.5,color='gainsboro')
        ax2.grid(True,which='both',lw=0.5,color='gainsboro')
        ax3.grid(True,which='both',lw=0.5,color='gainsboro')
        
        if logplot:
            for ax_ in [ax,ax2,ax3]:
                ax_.set_ylim([1e21,1e29])
            #ax2.ylim([1e20,9e29])
            #ax3.ylim([1e20,9e29])            
        else:
            for ax_ in [ax,ax2,ax3]:
                ax_.set_ylim([1e17,1e27])
               
        
        if plotMK == True:
            for ax_ in [ax,ax2,ax3]:
                ax_.set_xlim([min(data1['ts']),9])
            #ax2.xlim([min(data1['ts']),9])
            #ax3.xlim([min(data1['ts']),9])            
        else:
            for ax_ in [ax,ax2,ax3]:
                ax_.set_xlim([min(data1['ts']),max(data1['ts'])])
            #ax2.xlim([min(data1['ts']),max(data1['ts'])])
            #ax3.xlim([min(data1['ts']),max(data1['ts'])])
            
        ax.set_title('Compare DEMs - 1: '+title1+', 2: '+title2, fontsize = 30)
        ax2.set_title('1: '+title1, fontsize = 30)
        ax3.set_title('2: '+title2, fontsize = 30) 
        if peak_inset or high_inset:
            ax.legend(fontsize = 20, loc='upper left',markerscale=2)
        else:
            ax.legend(fontsize = 20, loc='lower left',markerscale=2)
            
        for ax_ in [ax2,ax3]:
            ax_.legend(ncol=2, fontsize = 20, loc='lower left',markerscale=2)
            
        if plotMK == True:
            for ax_ in [ax,ax2,ax3]:
                ax_.set_xlabel('T [MK]', fontsize=20)
                ax_.set_ylabel('$Emission\;Measure\;[cm^{-5}]$', fontsize=30)
                ax_.tick_params(axis='y', labelsize=25) #set_yticks(fontsize=25)
                ax_.tick_params(axis='x', labelsize=25)
        else:
            for ax_ in [ax,ax2,ax3]:
                ax_.set_xlabel('$\mathrm{\log_{10}T\;[K]}$', fontsize=30)
                ax_.set_ylabel('$Emission\;Measure\;[cm^{-5}]$', fontsize=30)
                ax_.tick_params(axis='y', labelsize=25) #set_yticks(fontsize=25)
                ax_.tick_params(axis='x', labelsize=25)
        
        #Get residuals, ready for comparison plot.
        residuals1, indices1, residuals2, indices2, inst = comparison_instruments(data1, data2)
        
        nf = len(inst)
        
        indices1 = [i-0.05 for i in indices1]
        indices2 = [i+0.05 for i in indices2]
        
        ax1.errorbar(indices1,residuals1, yerr=data1['edn'], fmt='^', color=fill1,
                     label='1: $\chi^2 =$ {:0.2f}'.format(chi_sq1) )#+' '+edn_string1)
        
        ax1.errorbar(indices2,residuals2, yerr=data2['edn'], fmt='^', color=fill2,
                     label='2: $\chi^2 =$ {:0.2f}'.format(chi_sq2) )#+' '+edn_string2)
        
        #Line at y=1 (ideal)
        ax1.plot([-1,nf+1],[1,1],'--',color='grey')

        ax1.set_ylim([0.5,2])
        ax1.set_xlim([-0.5,nf])
        ax1.set_yscale('log')
        ax1.set_title('DEM Residuals', fontsize=30)
        ax1.set_xticks(np.arange(nf),inst,rotation=45, fontsize=25, ha='right')
        ax1.tick_params(axis='x', length=20)
        ax1.set_yticks([0.5,1,2], ['0.5','1','2'], fontsize=25)
        #ax.set_xlabel('Channel', fontsize=30)
        ax1.set_ylabel('DN$_\mathrm{(DEM\ predicted)}$/DN$_\mathrm{(measured)}$', fontsize=30)
        ax1.yaxis.set_minor_formatter(NullFormatter())
        ax1.legend(fontsize = 20,markerscale=2)
        
        if tempcomp:
            plt.savefig('./'+timestring1+'/DEM_comparison_4panel_'+title1+'_vs_'+title2+'temp_range_comp.png')
        else:
            plt.savefig('DEM_comparison_4panel_'+title1+'_vs_'+title2+'.png')
            
#     m1 = DEMmax(data1['ts'], data1['DEM'], wind=3, plot=False)
#     m2 = DEMmax(data2['ts'], data2['DEM'], wind=3, plot=False)

#     print('DEM 1 is a maximum at: log(T)=', m1, 'OR, ', 10**m1/1e6, ' MK')
#     print('DEM 2 is a maximum at: log(T)=', m2, 'OR, ', 10**m2/1e6, ' MK')
        
    return consistent





def comparison_colors(data):
    """
    Makes a colors array in a way that will result in:
    -The same colors for each AIA channel, even if different numbers of channels are used
        (e.g. not the full 6).
    -Pull from the same pool of three for the colors for each XRT filter
    -Pull from the same pool of four for the colors for each NuSTAR energy range
    -Colors in the correct order.

    EXPECTS channels to be labeled as is output from dodem.
    """
    aia_colors = ['darkgreen', 'darkcyan', 'gold', 'sienna', 'indianred', 'darkorange']
    xrt_colors = ['darkslateblue', 'dodgerblue', 'cornflowerblue']
    nu_colors = ['purple', 'mediumpurple', 'plum', 'pink']
    
    expected_AIA = ['A94', 'A131', 'A171', 'A193', 'A211', 'A335']
    
    #If one of the channels is one of the expected aia channels, than add its color
    clrs=[]
    xrtcount=0
    nucount=0
    for i in range(0, len(data['chanax'])):
        chan = data['chanax'][i]
        if chan in expected_AIA:
            cc = [aia_colors[j] for j in range(0, len(expected_AIA)) if expected_AIA[j] == chan]
            clrs.extend(cc)
        if "med" in chan or "thin" in chan or "poly" in chan:
            clrs.append(xrt_colors[xrtcount])
            xrtcount+=1
        if "keV" in chan:
            clrs.append(nu_colors[nucount])
            nucount+=1
            
    return clrs     
    
def comparison_instruments(data1, data2):
    """
    Helps with the annoyances involved when comparing two DEMs done with different instrument 
    combinations.
    
    Make list of ALL involved instruments (between both results): AIA first, then XRT, then NuSTAR.
    Then, makes residuals arrays for each result so we can plot them together (no value included for
    cases where one result is missing an instrument that the other one does have.)
    
    """

    
    inst = []
    
    #AIA:
    expected_AIA = ['A94', 'A131', 'A171', 'A193', 'A211', 'A335']
    for i in range(0, len(data1['chanax'])):
        chan = data1['chanax'][i]
        if chan in expected_AIA:
            inst.append(chan)
    for i in range(0, len(data2['chanax'])):
        chan = data2['chanax'][i]
        if chan in expected_AIA and chan not in inst:
            inst.append(chan)
        
    #XRT
    for i in range(0, len(data1['chanax'])):
        chan = data1['chanax'][i]
        if "med" in chan or "thin" in chan or "poly" in chan:
            inst.append(chan)
    for i in range(0, len(data2['chanax'])):
        chan = data2['chanax'][i]
        if "med" in chan or "thin" in chan or "poly" in chan:
            if chan not in inst:
                inst.append(chan)
                
    #NuSTAR
    for i in range(0, len(data1['chanax'])):
        chan = data1['chanax'][i]
        if "keV" in chan:
            inst.append(chan)
    for i in range(0, len(data2['chanax'])):
        chan = data2['chanax'][i]
        if "keV" in chan:
            if chan not in inst:
                inst.append(chan)
    
    #We now have a list of all "instruments" that will be included in the residuals plot.
    #Later, we will want to plot the included instruments for each DEM, leaving blank spots 
    #when a given DEM did not use a particular instrument. For each DEM, make a list of
    #indices in the full mutual instrument list where that DEM did use a given instrument.
    
    indices = np.arange(len(inst))
    
    indices1 = [indices[i] for i in range(0, len(inst)) if inst[i] in data1['chanax']]
    indices2 = [indices[i] for i in range(0, len(inst)) if inst[i] in data2['chanax']]

    residuals1 = data1['dn_reg']/data1['dn_in']
    residuals2 = data2['dn_reg']/data2['dn_in']
    
    return residuals1, indices1, residuals2, indices2, inst
        
    


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)

def DEMmax(ts, DEM, wind=2, plot=False):
    """
    Does a lil Gaussian fit to find the maximum of the DEM.
    
    Keywords
    ---------
    
    ts, DEM - temperature and DEM arrays
    
    wind - Gaussian fit window: gives number of bins on each side of the 
            max value to use for the fit.
            
    """

    maxval=np.max(DEM)
    #Index of DEM maximum - first case if multiple temperature bins have same max value.
    pind = np.where(DEM == maxval)[0][0]
    maxtemp=ts[pind]
    dt = ts[pind]-ts[pind-1]
    fitrange = DEM[pind-wind:pind+wind+2]
    fitts = ts[pind-wind:pind+wind+2]
    try:
        popt, _ = optimize.curve_fit(gaussian, fitts, fitrange/maxval, p0=[1, maxtemp, dt])
    except ValueError:
        print('Max value possibly too close to the edge of the DEM array for Gaussian fit.')
        print('Max: ', maxval)
        print('DEM: ', DEM)
        print('Fitrange: ', fitrange)
        print('Just returning location of max value.')
        return (maxtemp, False)
        
    
    if plot:
        fig = plt.figure(figsize=(5, 5), tight_layout = {'pad': 1})
        plt.plot(ts, DEM/maxval)
        plt.plot(fitts,fitrange/maxval)
        plt.axvline(maxtemp, color='Red', linestyle='dashed', label='DEM Max')
        plt.axvline(popt[1], color='Red', label='Fit Peak')
        plt.plot(fitts, gaussian(fitts, *popt))
        plt.legend()
    
    return (popt[1], True)

    

def hightemp_EM(dem, ts, thresh, extract_vals=False, lowtemp_EM=False):
    """
    Take in temperature and dem arrays and calculate the approximate total emission predicted 
    above a certain threshold temperature (in log10(T)).
    
    """
    #Make a resampled temperature array with finer-spaced values
    newts = np.arange(ts[0], ts[-1], 0.01)
    newdem = np.interp(newts,ts,dem)
    
    if lowtemp_EM:
        lowT = np.where(newts < thresh)
        nTs = newts[lowT[0]]
        nDEM = newdem[lowT[0]]
    else:
        highT = np.where(newts > thresh)
        nTs = newts[highT[0]]
        nDEM = newdem[highT[0]]
    #print('Actual Threshold: logT =', round(nTs[0], 2), 'or', round(10**nTs[0]), 'K')
    
    
    res = np.trapz(nDEM, x=nTs)
    
    if extract_vals:
        return res, nTs
    else:
        return 'EM Above LogT='+str(round(nTs[0], 2))+f': {res:.2e}'+' cm^(-5)'

def both_powerlaws(ts, DEM, upper=True, lower=True, plot=True, fixlowerbound=False):
    """
    For given DEM solution, find upper and lower power law slopes of rise./decay around peak. 
    
    Lower Power Law Fit Bounds:
    ============================
    fixlowerbound=True: upper boundary is temp. closest to logT=6.35
    fixlowerbound=False: upper boundary is index of DEM max + 1
    
    """

    #Index of DEM maximum - first case if multiple temperature bins have same max value.
    maxdex = np.where(DEM == np.max(DEM))[0][0]
    
    
    if plot:
        fig = plt.figure(figsize=(5, 5), tight_layout = {'pad': 1})
        plt.semilogy(ts, DEM)
        plt.semilogy(ts[maxdex:], DEM[maxdex:])
    
    powerlaws = []

    if lower==True:
        
        #Find temp bin closest to 1 MK (minus 2 bins)
        from1mk = abs(ts-6)
        lowind = np.where(from1mk == np.min(from1mk))[0][0]-2
        
        
        from635 = abs(ts-6.35)
        upind = np.where(from635 == np.min(from635))[0][0]
        #print(ts[lowind])
        
        if fixlowerbound==False:
            upind=maxdex+1
        
        
        #Fit DEM from ~1MK to peak
        
        #Get into log-log space for linear fit
        ydata = np.log10(DEM[lowind:upind])
        xdata = ts[lowind:upind]
        res, cov = np.polyfit(xdata, ydata, 1, cov=True)
        m, b = res
        me = np.sqrt(np.diag(cov))[0]
        #print('Powerlaw: ', m)
        #print('Error: ', me)
        powerlaws.append((m, me, b))
        
        if plot:
            plt.semilogy(xdata, 10**(xdata*m+b))

    
    
    if upper==True:
        #Fit DEM from peak to max temp
        
        #Get into log-log space for linear fit
        ydata = np.log10(DEM[(maxdex+1):])
        xdata = ts[(maxdex+1):]
        res, cov = np.polyfit(xdata, ydata, 1, cov=True)
        m, b = res
        me = np.sqrt(np.diag(cov))[0]
        #print('Powerlaw: ', m)
        #print('Error: ', me)
        powerlaws.append((m, me, b))
        if plot:
            plt.semilogy(xdata, 10**(xdata*m+b))
    
        
    return powerlaws    
    
      
def get_DEM_params(file, save_params_file=False):
    
    """
    For a given dem results file, fetch assorted DEM result parameters 
    (written to work with output files named in the style used by dodem.high_temp_analysis(), where 
    the maximum temperature (maxT) is an optional input - default logT=7.2). 
    
    Keywords:
    -----------
    methodstring - 'iterative' - looks for outputs from xrt_dem_iterative2.pro wrapper
                '' - looks for outputs from DEMREG wrapper
                'MC' - looks for outputs from DEMREG wrapper (MCMC runs)
    
    
    Returns:
    -----------
    
    -m1: location of DEM peak (full temperature interval) 
    -max1: value of peak
    -above(i): DEM-estimated total emission measure above each temperature (in log(T)=x)
                i.e. above7 gives DEM-estimated EM above 10 MK.

    -DEM inputs (value, error, label) (dn_in, edn_in, chanax)

    -powerlaws: fit power law index+uncertainty both above and below the DEM peak
    
    """
    

    res = load_DEM(file)
    
    if res is not None:
        data, timestring, time = res
    else:
        return 
    
    
    #MAIN DEM
    lowdem = (data['DEM']-np.array(data['edem'])[0,:])
    hidem = (data['DEM']+np.array(data['edem'])[1,:])
    dem = data['DEM']
    ts = data['ts']
    
    m1, condition = DEMmax(data['ts'], dem, wind=3, plot=False)
    if condition==False:
        print(time)
    #print('DEM is a maximum at: log(T)=', m1, 'OR, ', 10**m1/1e6, ' MK')
    max1 = np.max(dem)
    
    above_peak, nTs = hightemp_EM(dem, ts, m1, extract_vals=True)
    below_peak, nTs = hightemp_EM(dem, ts, m1, extract_vals=True, lowtemp_EM=True)
    
    above_635, nTs = hightemp_EM(dem, ts, 6.35, extract_vals=True)
    below_635, nTs = hightemp_EM(dem, ts, 6.35, extract_vals=True, lowtemp_EM=True)

    
    
    above5, nTs = hightemp_EM(dem, ts, 6.7, extract_vals=True)
    above7, nTs = hightemp_EM(dem, ts, 6.84, extract_vals=True)
    above10, nTs = hightemp_EM(dem, ts, 7.0, extract_vals=True)
    
    above5l, nTs = hightemp_EM(lowdem, ts, 6.7, extract_vals=True)
    above7l, nTs = hightemp_EM(lowdem, ts, 6.84, extract_vals=True)
    above10l, nTs = hightemp_EM(lowdem, ts, 7.0, extract_vals=True)
    
    above5h, nTs = hightemp_EM(hidem, ts, 6.7, extract_vals=True)
    above7h, nTs = hightemp_EM(hidem, ts, 6.84, extract_vals=True)
    above10h, nTs = hightemp_EM(hidem, ts, 7.0, extract_vals=True)
    
    above5_ = [above5, above5l, above5h]
    above7_ = [above7, above7l, above7h]
    above10_ = [above10, above10l, above10h]
    
    
    powerlaws = both_powerlaws(ts, dem, plot=False, fixlowerbound=True)
    
    EMT_all = sum(dem*(10**ts))/sum(dem)/1e6
    index=14
    #print('Calculating EM-weighted T above:', data1['ts'][index])
    #thresh=5
    logthresh=6.7
    above_product, nTs = hightemp_EM(dem*(10**ts), ts, logthresh, extract_vals=True)    
    EMT_thresh= above_product/above5/1e6

    
    res = (m1, max1, above5_, above7_, above10_, 
           above_peak, below_peak, above_635, below_635,
           data['chanax'], data['dn_in'], data['edn_in'], powerlaws, EMT_all, EMT_thresh)

    save_params_file=True
    if save_params_file:

        savefile = file.split('.p')[-2]+'_withparams.pickle'
        #print(savefile)

        data['max'] = max1
        data['max_temp'] = m1
        data['above_5MK'] = above5_
        data['above_7MK'] = above7_
        data['above_10MK'] = above10_
        data['above_peak'] = above_peak
        data['below_peak'] = below_peak
        data['powerlaws'] = powerlaws
        data['EMT_all'] = EMT_all
        data['EMT_thresh_5'] = EMT_thresh

        with open(savefile, 'wb') as f:
             # Pickle the 'data' dictionary using the highest protocol available.
             pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)   
    
    return res      
    

    
def get_DEM_timeseries(time_intervals, working_dir, minT, maxT, name):
    
    """
    For a list of time intervals (and other characteristic parameters, allowing us to find your existing DEM output 
    pickle files), get the DEM output parameters for every file and compile into
    lists which are then saved to a nice output dictionary.
    
    """
    
    peaks=[]
    maxes=[]
    above5s=[]
    above7s=[]
    above10s=[]
    chanaxs=[]
    dn_ins=[]
    edn_ins=[]
    chanax_his=[]
    dn_in_his=[]
    edn_in_his=[]
    low_powers=[]
    hi_powers=[]
    above_peaks=[]
    below_peaks=[]
    above_635s=[]
    below_635s=[]    
    EMT_alls=[]
    EMT_threshs=[]

    result_time_intervals=[]


    for t in time_intervals:
        timestring=make_timestring(t)
        #print(timestring)
        file=working_dir+\
            timestring+'/'+timestring+'_'+str(minT)+'_'+str(maxT)+'_'+name+'_MC_DEM_result.pickle'
        params=get_DEM_params(file)
        if params is None:
            #Moving on if there's no DEM result file
            continue

        m1, max1, above5, above7, above10, \
           above_peak, below_peak, above_635, below_635,\
           chanax, dn_in, edn_in, powerlaws, EMT_all, EMT_thresh = params

        if len(dn_in) <= 6:
            #Not including cases of AIA-only DEMs when something went wrong with NuSTAR
            continue

        else:
            result_time_intervals.append(t)

        peaks.append(m1)
        maxes.append(max1)
        above5s.append(above5)
        above7s.append(above7)
        above10s.append(above10)
        chanaxs.append(chanax)
        dn_ins.append(dn_in)
        edn_ins.append(edn_in)
        low_powers.append(powerlaws[0][0])
        hi_powers.append(powerlaws[1][0])
        above_peaks.append(above_peak)
        below_peaks.append(below_peak)
        above_635s.append(above_635)
        below_635s.append(below_635)        
        EMT_alls.append(EMT_all)
        EMT_threshs.append(EMT_thresh)  


    vals = {'peaks': peaks,
        'maxes': maxes,
        'above5s': above5s,
        'above7s': above7s,
        'above10s': above10s,
        'chanaxs': chanaxs,
        'dn_ins': dn_ins,
        'edn_ins': edn_ins,
        'low_powers': low_powers,
        'hi_powers': hi_powers,
        'below_peaks': below_peaks,
        'above_peaks': above_peaks,
        'below_635s': below_635s,
        'above_635s': above_635s,                            
        'EMT_alls': EMT_alls,
        'EMT_threshs': EMT_threshs,
        'result_time_intervals': result_time_intervals}
    
    return vals

    
    
def pretty_orbit_timeseries(time_intervals, quantity, quantitylabel, label, color, backcolors,
                           error=False, quantity_low=[], quantity_high=[], errorcolor='gray',
                           comparisonbar=False, ylog=False, working_dir='./',
                           comp_band=[1.8e22, 1.5e23, 'Ishikawa (2017)'], plot_flares=False, 
                            show=True): 
    
    
    lw=2
    
    times = [t[0].datetime for t in time_intervals]
    midtimes=[(t[0]+(t[1]-t[0]).to(u.s)/2).datetime for t in time_intervals]


    starttime = (time_intervals[0][0]-120*u.s).datetime
    stoptime = (time_intervals[-1][1]+120*u.s).datetime

    times_ = copy.deepcopy(times)
    times_.append(time_intervals[-1][1].datetime)

    times2 = [t[1].datetime for t in time_intervals]
    
    fig, ax = plt.subplots(figsize=(15,4), tight_layout = {'pad': 1})

    flip=1
    for t in range(0, len(times)):
        if flip==1:
            ax.axvspan(times[t], times2[t], alpha=.5, color=backcolors[0])
        else:
            ax.axvspan(times[t], times2[t], alpha=.5, color=backcolors[1])
        flip*=-1

    if plot_flares:
        import pandas as pd
        import astropy.time
        df = pd.read_csv('fpmA.csv')
        starts = df['flare_start'].values
        stops = df['flare_end'].values
        

        early_starts = [(astropy.time.Time(s)-2*u.min).datetime for s in starts]
        late_stops = [(astropy.time.Time(s)+2*u.min).datetime for s in stops]  

        for j in range(0, len(early_starts)):
            ax.axvspan(early_starts[j], late_stops[j], alpha=.5, color='grey')
        
        
    if error:
        try:
            lp = np.hstack([quantity_low, quantity_low[-1]])
            hp = np.hstack([quantity_high, quantity_high[-1]])

            fill = ax.fill_between(times_, lp, hp, step="post", 
                                 color=color, alpha=0.1) 
        except ValueError:
            print('To plot error range, you need to define low/high bounds through the')
            print('quanity_low and quanity_high keywords.')
        
    ax.stairs(np.array(quantity), times_, linewidth=lw, color=color, 
              label=label,
             baseline=None)
    ax.set_ylabel(quantitylabel, fontsize=15, color=color)
    ax.tick_params(axis='y', labelcolor=color)
    ax.set_xlim([starttime, stoptime])
    ax.set_xlabel('Time', fontsize=15)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=1))
    
    if comparisonbar:
        ax.axhspan(comp_band[0], comp_band[1], color='brown', alpha=0.5, label=comp_band[2])   
        ax.legend()
    if error:
        if comparisonbar:
            themax = np.max([1.1*np.max(np.array(quantity_high)), 2*comp_band[1]])
            themin = np.min([0.5*np.min(np.array(quantity_low)), 0.5*comp_band[0]])
        else:
            themax = 1.1*np.max(np.array(quantity_high))
            themin = 0.5*np.min(np.array(quantity_low))
        
    else:
        if comparisonbar:
            themax = np.max([1.05*np.max(np.array(quantity)), 1.1*comp_band[1]])
            themin = np.min([0.95*np.min(np.array(quantity)), 0.5*comp_band[0]])    
        else:
            themax = 1.1*np.max(np.array(quantity))
            themin = 0.95*np.min(np.array(quantity))            

        
    ax.set_ylim([themin, themax])
    ax.set_title(label, fontsize=20, color=color)     
    
    if ylog:
        ax.set_yscale('log')
          
    savelabel= "-".join(label.split())        
        
    starttime = make_timestring(time_intervals[0])        
    plt.savefig(working_dir+'/'+savelabel+'_startat_'+starttime+'.png')

    if not show:
        plt.close()
    
    
    
    
def multi_orbit_summary(all_time_intervals, working_dir, name, minT=5.6, maxT=7.2):

    """
    Make a comparison plot for several DEM parameters over several orbits. 
    """
    
    #Plot inputs: edit here to add a new parameter
    valstrings = ['peaks','above10s','above7s','above5s', 'low_powers', 'hi_powers']
    backcolors = [['pink', 'lavenderblush'], ['powderblue', 'aliceblue'], ['powderblue', 'aliceblue'],
                 ['powderblue', 'aliceblue'], ['khaki', 'lemonchiffon'],['khaki', 'lemonchiffon']]
    colors = ['Red', 'Blue', 'Green', 'Purple', 'Orange', 'Red']
    factor = [1,1,1,1,1,-1]
    error = [False, True, True, True, False, False]
    theylabels=['T (MK)','EM (cm^-5)', 'EM (cm^-5)','EM (cm^-5)','Index', 'Index']
    thetitles=['DEM Peak Temperature', 'Total EM >10 MK', 'Total EM >7 MK', 'Total EM >5 MK',
              'Lower Power Law', 'Upper Power Law']    
    
 
    #list of dictionaries of DEM output values
    all_vals=[]
    #list of start of all time intervals for each orbit (datetime format)
    all_time_lists_nostep=[]
    #list of start of all time intervals + end of last one for each orbit (datetime format)
    all_time_lists = []
    #list of middle of all time intervals for each orbit (datetime format)
    all_midtimes = []
    #all start times (big list)
    alltimes = []
   
    #start, end for each orbit
    orbitstarts_ = []
    orbitstops_ =[]
   
    
    
    for ti in all_time_intervals:
       
        vals = get_DEM_timeseries(ti, working_dir, minT, maxT, name) 
        all_vals.append(vals)
        
        #different types of time list:
        #datetimes (start time of each interval)
        times = [t[0].datetime for t in ti]
        all_time_lists_nostep.append(times)
        alltimes.extend(times)
        #datetimes (same but with the end of the last interval added for plt.stairs use)
        times_ = copy.deepcopy(times)
        times_.append(ti[-1][1].datetime)
        all_time_lists.append(times_)
        #middle of each time interval
        midtimes=[(t[0]+(t[1]-t[0]).to(u.s)/2).datetime for t in ti]
        all_midtimes.append(midtimes)
       
        orbitstarts_.append(ti[0][0])
        orbitstops_.append(ti[-1][1])
       
    orbitstarts = [o.datetime for o in orbitstarts_]
    orbitstops = [o.datetime for o in orbitstops_]
 
 
 
    #MAKE BIG FIGURE
   
    lw=2.5
    tickfont=12
   
    fig, axes = plt.subplots(6, 1, figsize=(18,18), tight_layout = {'pad': 1}, sharex=True)
    fig.subplots_adjust(hspace=0.07)
   

    num=0
    while num < len(valstrings):
        ax=axes[num]
       
        for t in range(0, len(all_vals)):
            ax.axvspan(orbitstarts[t], orbitstops[t], alpha=.5, color=backcolors[num][0])
            if t > 0:
                ax.axvspan(orbitstops[t-1], orbitstarts[t], alpha=.5, color=backcolors[num][1])
       
        for t in range(0,len(all_vals)):
            plotvals = all_vals[t][valstrings[num]]
            plotvals = [factor[num]*v for v in plotvals]
           
            if error[num]:
                theval=np.array(plotvals)[:,0]
                quantity_low=np.array(plotvals)[:,1]
                quantity_high=np.array(plotvals)[:,2]
               
                ax.stairs(theval, all_time_lists[t], linewidth=lw, color=colors[num], baseline=None)
               
                try:               
                    lp = np.hstack([quantity_low, quantity_low[-1]])
                    hp = np.hstack([quantity_high, quantity_high[-1]])
   
                    
                    fill = ax.fill_between(all_time_lists[t], lp, hp, step="post",
                                    color=colors[num], alpha=0.1)
                except ValueError:
                    print('To plot error range, this needs to be a quantity with three values per timestep:')
                    print('(value, low, high)')
               
            else:
                ax.stairs(np.array(plotvals), all_time_lists[t], linewidth=lw, color=colors[num], baseline=None)
   
            if valstrings[num] == 'above10s':
                comp_band=[1.8e22, 1.5e23, 'Ishikawa et al. (2017)']
                ax.axhspan(comp_band[0], comp_band[1], color='brown', alpha=0.25, label=comp_band[2])
               
        ax.set_yscale('log')
        ax.set_title(thetitles[num], fontsize=20)
        ax.set_ylabel(theylabels[num], fontsize=15, color=colors[num])
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
   
        num+=1
   
    
    ax.xaxis.set_minor_locator(mdates.MinuteLocator(interval=10))
    ax.set_xlim([(np.min(orbitstarts_)-120*u.s).datetime, (np.max(orbitstops_)+120*u.s).datetime])
    ax.tick_params(axis='x', labelsize=tickfont)
    ax.tick_params(axis='y', labelsize=tickfont)
   
    
    wholestart=make_timestring([np.min(orbitstarts_).datetime, np.max(orbitstops_).datetime])
   
    plt.savefig(working_dir+'/'+wholestart+'_summary_plot.png')    
    
    
    
    
    
    
    
    