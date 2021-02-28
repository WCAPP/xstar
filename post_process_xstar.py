#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:07:30 2021

@author: patriciacho
"""

from astropy.table import Table
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib
import numpy as np
import os

matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['axes.linewidth'] = 0.5
matplotlib.rcParams['xtick.major.width'] = 0.4
matplotlib.rcParams['xtick.minor.width'] = 0.4
matplotlib.rcParams['ytick.major.width'] = 0.4
matplotlib.rcParams['ytick.minor.width'] = 0.4
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams.update({'font.size':10})

bohr = 5.292e-9      #bohr radius in cm
m_e=9.109e-28        #mass of electron
m_p=1.6726219e-24    #mass of proton
k_b=1.380648e-16
e=4.80e-10
planckev=4.135667516e-15
KtoeV = 11606
e_coulomb = 1.6021766e-19


# spectral_range = np.array([6.63,6.96])
# hnu_range = np.flip(12398.42/spectral_range,0)
# wide_range = np.array([4.3, 7.2])
# hnu_wide_range = np.flip(12398.42/wide_range,0)

# lxpos=0.95
# lypos=0.6

runs = ['run1_xstar2.56f','run1_xstar2.54a']
for run in runs:
    os.system('mkdir '+run+'/process_output')

#===========================================================================#
#                                                                           #
# Plot Charge State Distribution                                            #
#                                                                           #
#===========================================================================#

fig1 = plt.figure(figsize = [6,4])
ax1_0 = fig1.add_subplot(111)

for run in runs:
    abundances_fits = Table.read(run+'/xout_abund1.fits',hdu=1)
    
    radii = abundances_fits['radius']
    print(radii)
    print('number of zones: ', len(radii))
    print('log(xi): \n', abundances_fits['ion_parameter'])
    print('gas temperature in eV: \n', (abundances_fits['temperature']*1e4)/(KtoeV))
    
    elements = ['h_i','he_i','he_ii','li_i','li_ii','li_iii','be_i','be_ii','be_iii','be_iv','b_i','b_ii','b_iii','b_iv','b_v',
     'c_i','c_ii','c_iii','c_iv','c_v','c_vi','n_i','n_ii','n_iii','n_iv','n_v','n_vi','n_vii','o_i','o_ii','o_iii','o_iv',
     'o_v','o_vi','o_vii','o_viii','f_i','f_ii','f_iii','f_iv','f_v','f_vi','f_vii','f_viii','f_ix','ne_i','ne_ii','ne_iii',
     'ne_iv','ne_v','ne_vi','ne_vii','ne_viii','ne_ix','ne_x','na_i','na_ii','na_iii','na_iv','na_v','na_vi','na_vii',
     'na_viii','na_ix','na_x','na_xi','mg_i','mg_ii','mg_iii','mg_iv','mg_v','mg_vi','mg_vii','mg_viii','mg_ix','mg_x',
     'mg_xi','mg_xii','al_i','al_ii','al_iii','al_iv','al_v','al_vi','al_vii','al_viii','al_ix','al_x','al_xi','al_xii',
     'al_xiii','si_i','si_ii','si_iii','si_iv','si_v','si_vi','si_vii','si_viii','si_ix','si_x','si_xi','si_xii','si_xiii',
     'si_xiv','p_i','p_ii','p_iii','p_iv','p_v','p_vi','p_vii','p_viii','p_ix','p_x','p_xi','p_xii','p_xiii','p_xiv','p_xv','s_i',
     's_ii','s_iii','s_iv','s_v','s_vi','s_vii','s_viii','s_ix','s_x','s_xi','s_xii','s_xiii','s_xiv','s_xv','s_xvi','cl_i','cl_ii',
     'cl_iii','cl_iv','cl_v','cl_vi','cl_vii','cl_viii','cl_ix','cl_x','cl_xi','cl_xii','cl_xiii','cl_xiv','cl_xv','cl_xvi',
     'cl_xvii','ar_i','ar_ii','ar_iii','ar_iv','ar_v','ar_vi','ar_vii','ar_viii','ar_ix','ar_x','ar_xi','ar_xii','ar_xiii',
     'ar_xiv','ar_xv','ar_xvi','ar_xvii','ar_xviii','k_i','k_ii','k_iii','k_iv','k_v','k_vi','k_vii','k_viii','k_ix','k_x','k_xi',
     'k_xii','k_xiii','k_xiv','k_xv','k_xvi','k_xvii','k_xviii','k_xix','ca_i','ca_ii','ca_iii','ca_iv','ca_v','ca_vi','ca_vii',
     'ca_viii','ca_ix','ca_x','ca_xi','ca_xii','ca_xiii','ca_xiv','ca_xv','ca_xvi','ca_xvii','ca_xviii','ca_xix','ca_xx',
     'sc_i','sc_ii','sc_iii','sc_iv','sc_v','sc_vi','sc_vii','sc_viii','sc_ix','sc_x','sc_xi','sc_xii','sc_xiii','sc_xiv',
     'sc_xv','sc_xvi','sc_xvii','sc_xviii','sc_xix','sc_xx','sc_xxi','ti_i','ti_ii','ti_iii','ti_iv','ti_v','ti_vi','ti_vii',
     'ti_viii','ti_ix','ti_x','ti_xi','ti_xii','ti_xiii','ti_xiv','ti_xv','ti_xvi','ti_xvii','ti_xviii','ti_xix','ti_xx',
     'ti_xxi','ti_xxii','v_i','v_ii','v_iii','v_iv','v_v','v_vi','v_vii','v_viii','v_ix','v_x','v_xi','v_xii','v_xiii','v_xiv',
     'v_xv','v_xvi','v_xvii','v_xviii','v_xix','v_xx','v_xxi','v_xxii','v_xxiii','cr_i','cr_ii','cr_iii','cr_iv','cr_v','cr_vi',
     'cr_vii','cr_viii','cr_ix','cr_x','cr_xi','cr_xii','cr_xiii','cr_xiv','cr_xv','cr_xvi','cr_xvii','cr_xviii','cr_xix',
     'cr_xx','cr_xxi','cr_xxii','cr_xxiii','cr_xxiv','mn_i','mn_ii','mn_iii','mn_iv','mn_v','mn_vi','mn_vii','mn_viii','mn_ix',
     'mn_x','mn_xi','mn_xii','mn_xiii','mn_xiv','mn_xv','mn_xvi','mn_xvii','mn_xviii','mn_xix','mn_xx','mn_xxi','mn_xxii',
     'mn_xxiii','mn_xxiv','mn_xxv','fe_i','fe_ii','fe_iii','fe_iv','fe_v','fe_vi','fe_vii','fe_viii','fe_ix','fe_x','fe_xi',
     'fe_xii','fe_xiii','fe_xiv','fe_xv','fe_xvi','fe_xvii','fe_xviii','fe_xix','fe_xx','fe_xxi','fe_xxii','fe_xxiii',
     'fe_xxiv','fe_xxv','fe_xxvi','co_i','co_ii','co_iii','co_iv','co_v','co_vi','co_vii','co_viii','co_ix','co_x','co_xi',
     'co_xii','co_xiii','co_xiv','co_xv','co_xvi','co_xvii','co_xviii','co_xix','co_xx','co_xxi','co_xxii','co_xxiii',
     'co_xxiv','co_xxv','co_xxvi','co_xxvii','ni_i','ni_ii','ni_iii','ni_iv','ni_v','ni_vi','ni_vii','ni_viii','ni_ix','ni_x',
     'ni_xi','ni_xii','ni_xiii','ni_xiv','ni_xv','ni_xvi','ni_xvii','ni_xviii','ni_xix','ni_xx','ni_xxi','ni_xxii',
     'ni_xxiii','ni_xxiv','ni_xxv','ni_xxvi','ni_xxvii','ni_xxviii','cu_i','cu_ii','cu_iii','cu_iv','cu_v','cu_vi','cu_vii',
     'cu_viii','cu_ix','cu_x','cu_xi','cu_xii','cu_xiii','cu_xiv','cu_xv','cu_xvi','cu_xvii','cu_xviii','cu_xix','cu_xx',
     'cu_xxi','cu_xxii','cu_xxiii','cu_xxiv','cu_xxv','cu_xxvi','cu_xxvii','cu_xxviii','cu_xxix','zn_i','zn_ii','zn_iii',
     'zn_iv','zn_v','zn_vi','zn_vii','zn_viii','zn_ix','zn_x','zn_xi','zn_xii','zn_xiii','zn_xiv','zn_xv','zn_xvi','zn_xvii',
     'zn_xviii','zn_xix','zn_xx','zn_xxi','zn_xxii','zn_xxiii','zn_xxiv','zn_xxv','zn_xxvi','zn_xxvii','zn_xxviii','zn_xxix',
     'zn_xxx']
    
    print('ZONE ION FRACTIONS')
    
    ticker = 0
    for i in range(30):
        charges = np.zeros(0)
        fracs_beg = np.zeros(0)
        fracs_end = np.zeros(0)
        for j in range(i+1):
            info = abundances_fits[elements[ticker]]
            if (info[0]) != 0.0:
                print(info)
                charges = np.append(charges,j+1)
                fracs_beg = np.append(fracs_beg,float(info[0]))
                fracs_end=np.append(fracs_end,float(info[-1]))
            ticker += 1
        if len(charges) != 0:
            print('sum fractions: beginning, ending\n',sum(fracs_beg),sum(fracs_end))
            print('zbar: beginning, ending\n',sum(charges*fracs_beg),sum(charges*fracs_end))
            print('sigma_zbar: beginning, ending\n',sum(charges**2*fracs_beg)-
                  (sum(charges*fracs_beg))**2,sum(charges**2*fracs_end)-(sum(charges*fracs_end))**2)        
        if i==13:
            zbar_beg=sum(charges*fracs_beg)
            zbar_end=sum(charges*fracs_end)
    
            fig2 = plt.figure(figsize = [6,4])
            ax0 = fig2.add_subplot(111)
            ax0.plot(charges,fracs_beg,label=run+': beginning fractions')
            ax0.plot(charges,fracs_end,label=run+': ending fractions')
            ax0.set_title('Z = '+str(0.5*(zbar_beg+zbar_end))[:5]+
                          ', T = '+str(np.average((abundances_fits['temperature']*1e4)/(KtoeV)))[:5]+
                          r', log($\xi$) = '+str(np.average(abundances_fits['ion_parameter']))[:5]
                          )
            ax0.set_xlabel('charge')
            ax0.set_ylabel('fraction')
            ax0.xaxis.set_minor_locator(AutoMinorLocator())
            ax0.yaxis.set_minor_locator(AutoMinorLocator())
            ax0.tick_params(axis='both',which='both',direction='in')
            ax0.legend()
            plt.savefig(run+'/process_output/CSD.pdf')
            plt.close()
            
        
            ax1_0.plot(charges,fracs_beg,label=run+': beginning fractions')
            ax1_0.plot(charges,fracs_end,label=run+': ending fractions')
            
            
ax1_0.set_title('Z = '+str(0.5*(zbar_beg+zbar_end))[:5]+
              ', T = '+str(np.average((abundances_fits['temperature']*1e4)/(KtoeV)))[:5]+
              r', log($\xi$) = '+str(np.average(abundances_fits['ion_parameter']))[:5]
              )
ax1_0.set_xlabel('charge')
ax1_0.set_ylabel('fraction')
ax1_0.xaxis.set_minor_locator(AutoMinorLocator())
ax1_0.yaxis.set_minor_locator(AutoMinorLocator())
ax1_0.tick_params(axis='both',which='both',direction='in')
ax1_0.legend()
plt.savefig('CSD.pdf')
plt.close()    

#===========================================================================#
#                                                                           #
# Plot Spectra                                                              #
#                                                                           #
#===========================================================================#

fig1 = plt.figure(figsize = [5,7])
ax1_0 = fig1.add_subplot(211)
ax1_1 = fig1.add_subplot(212)
title = ''

for run in runs:
    emission_unit_conversion = 1./(0.36*0.6*2.*np.pi) * e_coulomb
    
    spect1_params = Table.read(run+'/xout_spect1.fits',hdu=1)
    spect1_xstar_spectra = Table.read(run+'/xout_spect1.fits',hdu=2)
                
    xstar_incident_spectrum = np.loadtxt('xs.dat',skiprows=1)
        
    x_new = spect1_xstar_spectra['energy']
    y_inward_new = spect1_xstar_spectra['emit_inward']
    y_outward_new = spect1_xstar_spectra['emit_outward']
    incident_spectrum = spect1_xstar_spectra['incident']
    transmitted_spectrum = spect1_xstar_spectra['transmitted']
    
    #===========================================================================#
    #                                                                           #
    # Check incident spectrum units are same as visrad                          #
    #                                                                           #
    #===========================================================================#
    
    fig = plt.figure(figsize = [6,4])
    ax0 = fig.add_subplot(111)
    ax0.plot(x_new,incident_spectrum*1e38*e_coulomb*(1./(4.*np.pi*(570.**2))),lw=3,label='xstar incident spectrum')
    ax0.plot(xstar_incident_spectrum[:,0],xstar_incident_spectrum[:,1],label='xs.dat',ls='--')
    ax0.legend()
    ax0.set_xlim(0,2500)
    ax0.set_xlabel('Energy [eV]')
    ax0.set_ylabel(r'Intensity $\left[\frac{J}{\rm{cm}^2\cdot \rm{s}\cdot \rm{eV}\cdot} \right]$')
    ax0.xaxis.set_minor_locator(AutoMinorLocator())
    ax0.yaxis.set_minor_locator(AutoMinorLocator())
    ax0.tick_params(axis='both',which='both',direction='in')
    plt.savefig(run+'/process_output/xstar_vs_visrad.pdf')
    plt.close()
    
    # x_old = t2_old['energy']
    # y_inward_old = t2_old['emit_inward']
    # y_outward_old = t2_old['emit_outward']
    
    fig = plt.figure(figsize=[5,7])
    ax0 = fig.add_subplot(211)
    ax0.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_inward_new,0),label='emitted inward', ls='--',color='blue',lw=0.75)
    ax0.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_outward_new,0),label='emitted outward',color='blue',lw=0.75)

    ax1_0.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_inward_new,0),label=run+': emitted inward', ls='--',lw=0.75)
    ax1_0.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_outward_new,0),label=run+': emitted outward',lw=0.75)
    
    # ax0.plot(x_old,y_inward_old,label='old version: emitted inward', ls='--',color='red',lw=0.75)
    # ax0.plot(x_old,y_outward_old,label='old version: emitted outward',color='red',lw=0.75)
    
    ax0.legend()
    ax0.set_xlim(1e-1,1e4)
    ax0.set_ylim(1e0,1e16)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylabel(r'Intensity $\left[\frac{J}{\rm{cm}^2\cdot \rm{s}\cdot \rm{eV}\cdot \rm{sr}}\right]$')
    ax0.tick_params(axis='both',which='both',direction='in')
    ax0.set_title('Z = '+str(0.5*(zbar_beg+zbar_end))[:5]+
                  ', T = '+str(np.average((abundances_fits['temperature']*1e4)/(KtoeV)))[:5]+
                  r', log($\xi$) = '+str(np.average(abundances_fits['ion_parameter']))[:5]
                  )    
    ax1 = fig.add_subplot(212)
    ax1.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_inward_new,0),label='emitted inward', ls='--',color='blue',lw=0.75)
    ax1.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_outward_new,0),label='emitted outward',color='blue',lw=0.75)
    
    # ax1.plot(x_old,y_inward_old,label='old version: emitted inward', ls='--',color='red',lw=0.75)
    # ax1.plot(x_old,y_outward_old,label='old version: emitted outward',color='red',lw=0.75)
    
    ax1.set_xlim(4.3,7.2)
    ax1.set_ylim(1e6,10**(12.8))
    ax1.set_yscale('log')
    ax1.set_ylabel(r'Intensity $\left[\frac{J}{\rm{cm}^2\cdot \rm{s}\cdot \rm{eV}\cdot \rm{sr}}\right]$')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(axis='both',which='both',direction='in')
    ax1.set_xlabel(r'Wavelength [$\AA$]')
    plt.savefig(run+'/process_output/emission_spectrum.pdf',bbox_inches='tight')
    plt.close()
    
    ax1_1.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_inward_new,0),label=run+': emitted inward', ls='--',lw=0.75)
    ax1_1.plot(12398.42/np.flip(x_new,0),1e38*emission_unit_conversion*np.flip(y_outward_new,0),label=run+': emitted outward',lw=0.75)
    title = title + 'Z = '+str(0.5*(zbar_beg+zbar_end))[:5]\
        +', T = '+str(np.average((abundances_fits['temperature']*1e4)/(KtoeV)))[:5]\
        +r', log($\xi$) = '+str(np.average(abundances_fits['ion_parameter']))[:5]+'\n'

ax1_0.legend()
ax1_0.set_xlim(1e-1,1e4)
ax1_0.set_ylim(1e0,1e16)
ax1_0.set_xscale('log')
ax1_0.set_yscale('log')
ax1_0.set_ylabel(r'Intensity $\left[\frac{J}{\rm{cm}^2\cdot \rm{s}\cdot \rm{eV}\cdot \rm{sr}}\right]$')
ax1_0.tick_params(axis='both',which='both',direction='in')
    
ax1_0.set_title(title[:-1])

ax1_1.set_xlim(4.3,7.2)
ax1_1.set_ylim(1e6,10**(12.8))
ax1_1.set_yscale('log')
ax1_1.set_ylabel(r'Intensity $\left[\frac{J}{\rm{cm}^2\cdot \rm{s}\cdot \rm{eV}\cdot \rm{sr}}\right]$')
ax1_1.xaxis.set_minor_locator(AutoMinorLocator())
ax1_1.tick_params(axis='both',which='both',direction='in')
ax1_1.set_xlabel(r'Wavelength [$\AA$]')
plt.savefig('emission_spectrum.pdf',bbox_inches='tight')
plt.close()  
