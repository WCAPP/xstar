pro xstar_read_output
; To learn about the extensions in a .fits file do:
;spect1_param=mrdfits(output_file, n,header)
;help, spect1_param
;print, header
;stop
;where n is from 0 to until the mrdfits returned value is 0
;pro xstar_read_output
; Duane uses 3ns for time integration, and 0.36*0.6cm2 for source area, and 4pi the later 2pi
emission_unit_conversion = 1./(0.36*0.6*2.*!pi) * !pcon.e ; such that xstar erg/s/erg multiplied by that factor gives J/cm2/eV/sr/s

xstar_version='_11_2018'
;xstar_version='231'

;CD, '/Users/gploise/Documents/Experimental/Sandia/ZAPP/Si_photo/astro/xstar/runs/run'+trim(run)+'/'
CD, '/Users/gploise/Documents/Experimental/Sandia/ZAPP/Si_photo/astro/xstar/runs_xstar'+xstar_version+'/'

runs = 40 + indgen(9)
;runs=7
runs=11
runs=20
runs=48
runs='_stephanie'
runs='_nlte10'
;runs=14
;runs=16
;runs='_kallman_Tc59eV_Tb56eV'

buffer=1

nruns=n_elements(runs)


visrad_times_with_extrapol=[-3.00, -2.02, -1.22, -0.42, 0.38, 1.18, 1.98, 2.78, 3.00] ; ns
ntimes=n_elements(visrad_times_with_extrapol)
times=(visrad_times_with_extrapol-visrad_times_with_extrapol[0])*1.e-9 ; seconds
delta_time = [times[1:ntimes-1], times[ntimes-1]]-[times[0],times[0:ntimes-2]] ;for time integration, trapeze scheme = sum_i fi*dti
integration_duration=last(times)-times[0]
weights=0.5*delta_time

if nruns ne 1 and nruns ne n_elements(weights) then stop
if nruns eq 1 then integration_duration = 3.e-9

; Z data wda
wdf_z_tra=aw(1)
file='/Users/gploise/Documents/Experimental/Sandia/ZAPP/Si_photo/Chris_Fontes/6_tra_compile_ave.txt'
data=read_ascii(file, COMMENT_SYMBOL='#', data_start=0)
r=execute(trim('data=data.'+tag_names(data)))
x=reform(data[0,*])
y=reform(data[1,*])
i2w, x, y, wdf_z_tra, /xy, cl='EXPT_tra_compile_ave', xl='Wavelength [Å]', yl='Transmission'
abs_resolving_power = 2300


;z_emission_data_file = '/Applications/Prism/SPECT3D_14.0.0/silicon/z2674_6mm_length.txt'
;emi_resolving_power=2600. ; spectro sept 2013
z_emission_data_file = '/Applications/Prism/SPECT3D_14.0.0/silicon/z2786_2787_ave_3mm_length.txt'
emi_resolving_power=4144. ; spectro april 2014
;z_emission_data_file = '/Applications/Prism/SPECT3D_14.0.0/silicon/z2672_12mm_length.txt'
;emi_resolving_power=2600.; spectro sept 2013
wdf_z_emi = aw(1)
data=read_ascii(z_emission_data_file, COMMENT_SYMBOL='#', data_start=0)
r=execute(trim('data=data.'+tag_names(data)))
x=reform(data[0,*])
y=reform(data[1,*])
i2w, x, y, wdf_z_emi, /xy, cl='Z data', xl='Wavelength [Å]', yl='Emergent intensity (a.u.)'
wdf_z_emi_rescale = aw(1)

wda_inc = aw(nruns)
wda_trans = aw(nruns)
wda_emi = aw(nruns)
wda_emi_time_int = aw(1)
wda_trans_time_int = aw(1)
wda_inc_time_int = aw(1)
wda_transmission_time_int = aw(1)

spectral_range=[6.63, 6.96]
hnu_range = 12398.42/spectral_range 
hnu_range = reverse(hnu_range)
wide_range = [4.3, 7.2]
hnu_wide_range = 12398.42/wide_range
hnu_wide_range = reverse(hnu_wide_range)


lxpos=0.95
lypos=0.6

for irun=0, nruns-1 do begin
    run=runs[irun]
    wda=aw(14)

    CD, 'run'+trim(run)+'/'
    
    output_file='xout_abund1.fits'
    abundances=mrdfits(output_file, 'ABUNDANCES',header)
    ; contains nuc charge / ion charge fields
    ; .ion_parameter =20.21
    ; .pressure
    ; .radius (1st radius)
    ; .temperature=66.03
    ; .x_e=1.02
    ;print,header
    print, 'ZONE PARAMETERS'
    radii=abundances.radius
    print,'radius in cm',radii
    print, '# of zones', n_elements(radii)
    print,'log(xi [erg.cm/s])',abundances.ion_parameter
    ;print,'delta radius in cm',radii-radii[0], format='(E12.3)'
    print,'gas temperature in eV',abundances.temperature*1.e4/(!pcon.e/!pcon.kb)
    ;print,'gas temperature in K',abundances.temperature*1.e4
    
    print, 'ZONE ION FRACTIONS'   
    ; print all non zero abundances
    for z=1, 30 do begin
      charges=[]
      fracs_first=[]
      fracs_last=[]
      ;  sum_frac=0.
      for charge=0, z-1 do begin
        zion2symb,z,charge+1,symb,ziform='Z_ION'
        symb=symb[0]
        ;    print,symb
        r=execute('frac=abundances.'+symb) 
        if frac[0] ne 0. then begin
          print,symb,frac
          charges=[charges,charge]
          fracs_first=[fracs_first,frac[0]]
          fracs_last=[fracs_last,last(frac)]
        endif
      endfor
      ;  print,charges,fracs
      ; save
      if charges ne !NULL then begin
        print,'sum fractions, first, last',total(fracs_first),total(fracs_last)
        print,'zbar, first, last',total(charges*fracs_first),total(charges*fracs_last)
        print,'sigma_zbar, first, last',total(charges^2*fracs_first)-(total(charges*fracs_first))^2,total(charges^2*fracs_last)-(total(charges*fracs_last))^2
        if strcmp(strmid(symb,0,2),'Si') then begin
          zbar_first=total(charges*fracs_first)
          zbar_last=total(charges*fracs_last)
          i2w,charges,fracs_first,wda[8],/xy, xl='ion', yl='fraction',cl='First, Z='+trim(zbar_first, '(f10.2)')
          i2w,charges,fracs_last,wda[9],/xy, xl='ion', yl='fraction',cl='Last, Z='+trim(zbar_last, '(f10.2)')     
          wriasc,strmid(symb,0,2)+'_csd_first.txt',wda[8],formatout=0
          wriasc,strmid(symb,0,2)+'_csd_last.txt',wda[9],formatout=0
          print,'*******'
          print,strmid(symb,0,2)+'_csd_first.txt', ' file written'
          print,strmid(symb,0,2)+'_csd_last.txt', ' file written'
        endif
      endif
    endfor
    
    
    
    csd_label = 'Z='+trim(0.5*(zbar_first+zbar_last), '(f10.2)')+', T='+trim(mean(abundances.temperature*1.e4/(!pcon.e/!pcon.kb)),'(f10.2)')+' eV, log($\xi$) = ' +trim(mean(abundances.ion_parameter),'(f10.2)')

    ; read spectral output
    print, 'READ SPECTRAL OUTPUT'
    output_file='xout_spect1.fits'
    
    spect1_param=mrdfits(output_file, 'PARAMETERS',header)
    spect1_xstar_spectra=mrdfits(output_file, 'XSTAR_SPECTRA',header)
    
    n = n_elements(spect1_xstar_spectra)
    hnu=dblarr(n)
    inc = hnu
    trans = hnu
    emit_in = hnu
    emit_out = hnu
    scatt = hnu
    for i=0L, n-1 do begin
      struct = spect1_xstar_spectra[i]
      hnu[i] = struct.energy
      inc[i] = struct.INCIDENT
      trans[i] = struct.TRANSMITTED
      emit_in[i] = struct.emit_inward
      emit_out[i] = struct.emit_outward
      ;  if i eq 2 then print, emit_out[i]
;      scatt[i] = struct.scattered
    endfor
    
    ; calculate spectral grid resolution in ROI
    ind=where((hnu ge hnu_range[0]) and (hnu le hnu_range[1]))
    hnu_res=mean(diff(hnu[ind]))
    lambda_res = abs(mean(diff(12398419./hnu[ind]))) ; mA
    print, 'hnu_res eV', hnu_res
    print, 'lambda_res mÅ', lambda_res
    
    ; put all spectra in wdf arrays
    index=where((hnu ge hnu_range[0]) and (hnu le hnu_range[1]))
    i2w, hnu[index], inc[index]*1.e38*emission_unit_conversion, wda[0], /xy, cl='incident', xl='Photon energy [eV]', yl='Intensity [J/cm$^2$/eV/sr/s]'
    i2w, hnu[index], trans[index]*1.e38*emission_unit_conversion, wda[1], /xy, cl='transmitted', xl='Photon energy [eV]', yl='Intensity [J/cm$^2$/eV/sr/s]'
    i2w, hnu[index], emit_in[index]*1.e38*emission_unit_conversion, wda[2], /xy, cl='emit_in', xl='Photon energy [eV]', yl='Intensity [J/cm$^2$/eV/sr/s]'
    i2w, hnu[index], emit_out[index]*1.e38*emission_unit_conversion, wda[3], /xy, cl='emit_out energy', xl='Photon energy [eV]', yl='Intensity [J/cm$^2$/eV/sr/s]'
    i2w, 12398.42/hnu[index], emit_out[index]*1.e38*emission_unit_conversion, wda[4], /xy, cl='emit_out wavelength', xl='Wavelength [Å]', yl='Intensity [J/cm$^2$/eV/sr/s]'
    i2w, hnu[index], scatt[index]*1.e38*emission_unit_conversion, wda[5], /xy, cl='scattered', xl='Photon energy [eV]', yl='Intensity [J/cm$^2$/eV/sr/s]'
    ;transmission
    i2w, hnu[index], trans[index]/inc[index], wda[6], /xy, cl='transmission energy', xl='Photon energy [eV]', yl='Transmission'
    i2w, 12398.42/hnu[index], trans[index]/inc[index], wda[7], /xy, cl='transmission wavelength', xl='Wavelength [Å]', yl='Transmission'

    index=where((hnu ge hnu_wide_range[0]) and (hnu le hnu_wide_range[1]))
    i2w, 12398.42/hnu[index], emit_out[index]*1.e38*emission_unit_conversion, wda[12], /xy, cl='emit_out wavelength wide', xl='Wavelength [Å]', yl='Intensity [J/cm$^2$/eV/sr/s]'

    if nruns eq 1 then mul, wda[4], 0, integration_duration , yl='Intensity [J/cm$^2$/eV/sr]' , cl='emit_out wavelength time int'
    
    ; apply instrumental convolution
    wdgaussconv, wda[7], wda[10], rp=abs_resolving_power
    wdgaussconv, wda[4], wda[11], rp=emi_resolving_power
    wdgaussconv, wda[12], wda[13], rp=emi_resolving_power
    
    ; check parameters or the incident spectrum in same unit as in the visrad input
    if file_test('xs.dat') then begin
      flux=inc *!pcon.e* 1.e38 ; J/s/eV
      flux=flux * 1./(4.*!pi*570.^2) ; J/s/eV/cm2
      data=read_ascii('xs.dat',data_start=1,delimiter=' ', num_records=200)
      x=reform(data.field1[0,*]); ev
      y=reform(data.field1[1,*]); J/cm2/s/eV
      plot, x, y, xra=[0., 2500.], title='run'+trim(run)
      oplot, hnu, flux, col=1, line=2
      spng, 'visrad_vs_xstar.png',res=100
    endif
    ; CSD
    p=plotw(wda[8:9], lypos=400./900., title=csd_label, lxpos=0.65, font_size=50, xtitle='', thick=[7,7],/buffer)
    p.save, 'CSD.png', res=100
    
    ; TRANSMISSION
    p=plotw([wdf_z_tra,wda[7]], xr=spectral_range,yrange=[0., 1.], title = 'run'+trim(run)+', xstar res='+trim(lambda_res,'(f10.3)')+' mÅ='+trim(hnu_res,'(f10.3)')+' eV', lxpos=lxpos, lypos=lypos, col=color([1,3]),buffer=buffer)
    p=plotw(wda[8:9], title=csd_label, lxpos=0.83,lypos=0.37, font_size=15, xtitle='', thick=[7,7], current=1, position=[0.7,0.2, 0.9,0.4], xrange=[7.,14.])
    p.save, 'transmission.png', res=100
    p=plotw([wdf_z_tra,wda[10]], xr=spectral_range,yrange=[0., 1.], title = 'run'+trim(run)+', xstar res='+trim(lambda_res,'(f10.3)')+' mÅ='+trim(hnu_res,'(f10.3)')+' eV, convolv='+trim(abs_resolving_power), lxpos=lxpos, lypos=lypos, col=color([1,3]),buffer=buffer)
    p=plotw(wda[8:9], title=csd_label, lxpos=0.83,lypos=0.37, font_size=15, xtitle='', thick=[7,7], current=1, position=[0.7,0.2, 0.9,0.4], xrange=[7.,14.])
    p.save, 'transmission_conv.png', res=100

    ; EMISSION
    p=plotw(wda[[0,1,2,3,5]], xr=hnu_range, title = 'run'+trim(run), thick=4+intarr(5), lxpos=lxpos, lypos=lypos,buffer=buffer);, yrange=[0., 1.e-2])); line=indgen(5))
    p.save, 'emission_all.png', res=100

    rescale_factor = mx(wda[11], left=6.70, right=6.75)/mx(wdf_z_emi, left=6.70, right=6.75)
    xfr, wdf_z_emi, wdf_z_emi_rescale
    mul, wdf_z_emi_rescale, 0, rescale_factor, cl=wdcomment(wdf_z_emi)+', rescale = '+trim(rescale_factor)   
    p=plotw([wdf_z_emi_rescale,wda[11]], xr=spectral_range, thick=[4,4], col=color([1,3]), lxpos=lxpos, lypos=lypos, title = 'run'+trim(run),buffer=buffer);, yrange=[0., 1.e-2])); line=indgen(5))
;    p=plotw([wdf_z_emi_rescale,wda[11], wda[4]], xr=spectral_range, thick=[4,4,2], col=color([1,3,2]), lxpos=lxpos, lypos=lypos, title = 'run'+trim(run),buffer=buffer);, yrange=[0., 1.e-2])); line=indgen(5))
    p.save, 'emission_out.png', res=100
        

    p=plotw(wda[12:13], xr=wide_range,/ylog, thick=[4,4], col=color([2,3]), lxpos=lxpos, lypos=lypos, title = 'run'+trim(run),buffer=buffer);, yrange=[0., 1.e-2])); line=indgen(5))
    p.save, 'emission_out2.png', res=100
  
 
    pff_file='xstar_run'+trim(run)+'.pff'
    createpff, pff_file, /delete
    for i=0, n_elements(wda)-1 do wri, wda[i], /full_precision
    dirp
    closepff
    
    if nruns ne 1 then begin
      xfr, wda[0], wda_inc[irun] 
      xfr, wda[1], wda_trans[irun] 
      xfr, wda[11], wda_emi[irun]
      wdgaussconv, wda_inc[irun], wda_inc[irun], rp=abs_resolving_power
      wdgaussconv, wda_trans[irun], wda_trans[irun], rp=abs_resolving_power
    endif
    
;    delwdf, wda
    cd, '..'
endfor
if nruns ne 1 then begin 
  weighted_sum, wda_inc, wda_inc_time_int , weights=weights
  weighted_sum, wda_trans, wda_trans_time_int , weights=weights
  div, wda_trans_time_int, wda_inc_time_int , wda_transmission_time_int
  putx, 12398.42/getxf(wda_transmission_time_int), wda_transmission_time_int


  weighted_sum, wda_emi, wda_emi_time_int , weights=weights
  rescale_factor = mx(wda_emi_time_int, left=6.70, right=6.75)/mx(wdf_z_emi, left=6.70, right=6.75)
  xfr, wdf_z_emi, wdf_z_emi_rescale
  mul, wdf_z_emi_rescale, 0, rescale_factor, cl=wdcomment(wdf_z_emi)+', rescale = '+trim(rescale_factor)

  p=plotw([wdf_z_tra, wda_transmission_time_int], xr=spectral_range,yrange=[0., 1.], col=color([1,3]), $
    title = 'xstar res='+trim(lambda_res,'(f10.3)')+' mÅ='+trim(hnu_res,'(f10.3)')+' eV, convolv='+trim(abs_resolving_power), lxpos=lxpos, lypos=lypos);, col=color([1,3]))
  p.save, 'run'+trim(runs[0])+'_time_int_tra.png', res=100
  p=plotw([wdf_z_emi_rescale, wda_emi_time_int], xr=spectral_range, ytitle='Intensity [J/cm$^2$/eV/sr]', col=color([1,3]), $
    title = 'xstar res='+trim(lambda_res,'(f10.3)')+' mÅ='+trim(hnu_res,'(f10.3)')+' eV, convolv='+trim(emi_resolving_power), lxpos=lxpos, lypos=lypos);, col=color([1,3]))
  p.save, 'run'+trim(runs[0])+'_time_int_emi.png', res=100
endif
end
