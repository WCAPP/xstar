;CD, '/Users/gploise/Documents/Experimental/Sandia/ZAPP/Si_photo/astro/xstar/runs_xstar221/'
CD, '/Users/gploise/Documents/Experimental/Sandia/ZAPP/Si_photo/astro/xstar/runs_xstar_11_2018/'

run=3


; create run specific folder
runfolder='run'+trim(run)+'/'
if ~file_test(runfolder) then  FILE_MKDIR, runfolder


; create spectrum file
file='/Users/gploise/Documents/Experimental/Sandia/ZAPP/Si_photo/astro/run57_ave_visrad_all_extrap.dat' ; time average of side pinch radaition project before 2016
;file='/Applications/Prism/PrismSPECT_silicon/PS_vs_XSTAR/xs_nlte10_3p.dat'; this one is one of the drives proposed for the NLTE-10 workshop
data=read_ascii(file,data_start=9,delimiter=' ', num_records=200)
x=reform(data.field1[0,*])
y=reform(data.field1[1,*])

hnu=genarray(5.,5000.,200.)
linterp, x, y, hnu, spect, missing='NaN'
N=n_elements(hnu)
openw, lun, runfolder+'xs.dat', /get_lun
printf, lun, N, format='(i4)'
for i=0,N-1 do begin
  printf, lun, hnu[i],spect[i], format='(e10.4," ", e9.3)'
endfor
close, lun
free_lun, lun

flux=int_tabulated(hnu, spect ,/sort, /double)*1.e7 ;= 1.0531e12*1.E7;erg/s/cm2 ; from visrad optimiz, time average
print, 'time average spectral radiance erg/s/cm2', flux

; my definitions out of experimental data, best estimates as of end of 2012
;density_si=8.25e17;at/cc
;mean_charge=10.4 ; from PS
;nne=density_si*mean_charge
;density=1.e9 ; hydrogen nucleus density n , if fully ionized and solar abund then particle density is 2.3n 
;
;flux=1.3e12*1.E7;erg/s/cm2 ; from visrad; gives a 59eV brightness temperature
;xi=4.*!pi*flux/(nne); phot param 
;trad=59. ; in eV , from flux
;temp=20.; eV, from line ratio
;radius=4.4; cm distance source to plasma
;slab_length=1.0;cm

; redefine parameters in the similar than Roberto with Cloudy input file


; #### SECTION defining all parameters I want to impose
temp=33.;eV
temp=30.;eV NLTE-0 workshop


; measured quantities
si_areal_density = 3.06e17
o_areal_density = 1.5e17
 
; plasma dimensions: 
; to calculate the Si density
plasma_size_fwhm_in_abs_los = 0.36; cm
plasma_size_fwhm_in_abs_los = 0.1; cm, NLTE-10 workshop

;plasma_size_fwhm_in_abs_los = 0.36/2. ; to increase density by a factor of 2
;plasma_size_fwhm_in_abs_los = 0.36/4. ; to increase density by a factor of 4


; Plasma size in the column density, should be plasma_size_fwhm_in_abs_los of rabsoprtion but may 
; be varied for simulating emission, ?ask Kallman.  
; before run 12:
;plasma_size_fwhm = 0.36 ; cm
plasma_size_fwhm = 0.3 ; cm ; run12
plasma_size_fwhm = 0.6 ; cm ; run13
plasma_size_fwhm = 1.2 ; cm ; run14
plasma_size_fwhm = plasma_size_fwhm_in_abs_los ; cm this si the size to use for simulating absorption

;density_si=1.07e18;at/cc
; The density spatial distribution is modeled as a top hat over the FWHM of the observed distribution. 
density_si = si_areal_density/plasma_size_fwhm_in_abs_los ;at/cc
density_si = 8.5e17 ;at/cc
print, 'density_si', density_si
;stop
radius_in = 570.;cm
radius_out = radius_in + plasma_size_fwhm;cm
slab_length = radius_out-radius_in

frac_Si_ovr_H = 10.

frac_O_ovr_Si = o_areal_density / si_areal_density
frac_O_ovr_H = frac_O_ovr_Si * frac_Si_ovr_H

density=density_si/frac_Si_ovr_H ; hydrogen nucleus density n , if fully ionized and solar abund then particle density is 2.3n 
print,'Hydrogen density',density


; define simulations variables
modelname='run'+trim(run)

cfrac=0.
temperature=temp*!pcon.k_11604*1./(1.e4); in 1e4 K ; if niter=0 then temp is fixed 
print, 'temperature in 1e4K', temperature
pressure=4.4e7 ;dyne/cm2 (1 dyne/cm2 = 0.1Pa) ignored if lcpres=0
lcpres=0 ; 0 cst density case and pressure ignored (and xi used), if 1 cst pressure (and Xi used)

; if we preserves , flux on sample, inner radius, and exp photoionization parameter (4 pi F/ ne), then  density =ne_exp which is far too large
; one can not preserves that many parameters
;print,'hydrogen density=',density

spectrum='"file"'
;trad=trad*!pcon.k_11604*1./(1.e7) ; in unit of 1.e7 K for bbody, or keV for Brem

rlrad38=flux*4.*!pi*radius_in^2./(1.E38) ;in 1.e38 erg/s,  luminosity integrated between 1 and 1000 Ry. 1.3E12 W/cm2 over 1by0.4 cm2 sample
print,'rlrad38=',rlrad38

rlogxi=alog10(4.*!pi*flux/density)  ; xi=L/(n.R^2) used to infer the inner most edge by XSTAR
print,'rlogxi=',rlogxi

if strcmp(spectrum,'"file"') then trad=0

column=density*slab_length  ; with the density, used to compute the slab thickness
print, 'H column', column


  ; elt frac abundances relaitove to solar (fraction between 0. (0%) and 100.(10 000%)), grevesse & sauval 1996
  habund=1.0
  heabund=1.0; Nov 2018: Tim told me that it needs to have Helium
  liabund=0.0
  beabund=0.
  babund=0.
  cabund=0.
  nabund=0.
  solar_abund_o=1./10.^(12.-8.83251)
  oabund=0.;frac_o_ovr_H/solar_abund_o
  print, 'oabund',oabund
  ;oabund=0.
  fabund=0.
  neabund=0.
  naabund=0.
  mgabund=0.
  alabund=0.
  solar_abund_si=1./10.^(12.-7.54407)
  siabund=frac_Si_ovr_H/solar_abund_si
  print, 'siabund',siabund
  pabund=0.
  sabund=0.
  clabund=0.
  arabund=0.
  kabund=0.
  caabund=0.
  scabund=0.
  tiabund=0.
  vabund=0.
  crabund=0.
  mnabund=0.
  feabund=0.
  coabund=0.
  niabund=0.
  cuabund=0.
  znabund=0.

  ;hidden variables
  niter=-99 ; for thermal equilibrium and charge conservation
  npass=3 ; must be odd and better if greater than 1
  critf=1.E-08
  nsteps=1 ; limit max spatial zone size
  xeemin=0.1
  emult=0.1
  taumax=20.
  lwrite=1
  lprint=1
  lstep=0
  ncn2=999 ; spectral resolution, higher means more spectral resoltuino, range= 999 to 999999, 99999 is decent for the transmission measurements
  radexp=0
  vturbi=0.
;  loopcontrol=0
;  mode='"ql"'
;  
  format_exponential = '(e15.4)'

  ; generate xstar case, write the comand line file
  file='./run'+trim(run)+'/run'+trim(run)+'.sh'
  openw,u,file,/get
  printf,u,'rm *.fits'
  printf,u,'xstar \'
  printf,u,'cfrac='+trim(cfrac,'(f5.3)')+' \'
  printf,u,'temperature='+trim(temperature,format_exponential)+' \'
  printf,u,'pressure='+trim(pressure)+' \'
  printf,u,'density='+trim(density,format_exponential)+' \'
  printf,u,'spectrum='+spectrum+' \'
  printf,u,'spectrum_file=xs.dat \'
  printf,u,'spectun=0 \'
  printf,u,'trad='+trim(trad,'(b1)')+' \'
  printf,u,'rlrad38='+trim(rlrad38,format_exponential)+' \'
  printf,u,'column='+trim(column,format_exponential)+' \'
  printf,u,'rlogxi='+trim(rlogxi,format_exponential)+' \'
  printf,u,'lcpres='+trim(lcpres,'(b1)')+' \'
  for z=1,30 do begin
    zion2symb,z,0,elt,ziform='Z' ; deifine elt name for a given nucleus charge
    r=execute('xabund='+elt[0]+'abund')
    printf,u,elt+'abund='+trim(xabund,'(e10.3)')+' \'
  endfor
  printf,u,'modelname='+modelname+' \'
  printf,u,'niter='+trim(niter)+' \'
  printf,u,'npass='+trim(npass)+' \'
  printf,u,'critf='+trim(critf)+' \'
  printf,u,'nsteps='+trim(nsteps)+' \'
  printf,u,'xeemin='+trim(xeemin)+' \'
  printf,u,'emult='+trim(emult)+' \'
  printf,u,'taumax='+trim(taumax)+' \'
  printf,u,'lwrite='+trim(lwrite,'(b1)')+' \'
  printf,u,'lprint='+trim(lprint,'(b1)')+' \'
  printf,u,'lstep='+trim(lstep,'(b1)')+' \'
  printf,u,'ncn2='+trim(ncn2)+' \'
  printf,u,'radexp='+trim(radexp)+' \'
  printf,u,'vturbi='+trim(vturbi)
  ;printf,u,'loopcontrol='+trim(loopcontrol)+' \'
  ;printf,u,'mode='+mode+' \'
  close,u
  free_lun,u
  spawn,'/usr/local/bin/bbedit '+file
  spawn,'chmod u+x '+file

end