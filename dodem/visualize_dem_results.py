import pickle
from matplotlib import pyplot as plt
import numpy as np
import copy

import matplotlib.dates as mdates
from matplotlib.ticker import NullFormatter, ScalarFormatter
from astropy import units as u
import os


def load_DEM(time, filename):
    """
    Take in a DEM file (saved output from dodem.dodem) + return dictionary of results + inputs.
    
    CURRENT FUNCTIONALITY: Enter name of file (pickle) with DEM results.
    
    OLD FUNCTIONALITY: Instead of file name, requests time interval and fpm (standardized file 
                        names+locations assumed).
    """
    
    timestring = time[0].strftime('%H-%M-%S')
    stopstring = time[1].strftime('%H-%M-%S')
    #print(timestring, stopstring)
    timestring=timestring+'_'+stopstring 
    
    file=filename
    
    with open(file, 'rb') as f:
        data = pickle.load(f)
        
    return data, timestring

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
    
    plt.close(fig)
    plt.savefig(title+'_DEM_plot.png')
    
    

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
            plt.close(fig)
            plt.savefig('./'+timestring1+'/DEM_comparison_'+title1+'_vs_'+title2+'temp_range_comp.png')
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
        
    
    
    
    
    
    
    
    
    