# -*- coding: utf-8 -*-

import os
import sys

input = sys.argv[1]

directory_name = str(input)

with open ('run_xstar.sh', 'w') as rsh:
    rsh.write('xstar \
        cfrac=0.000 \
        temperature=3.8295e+001 \
        pressure=4.4E+007 \
        density=8.5000e+016 \
        spectrum="file" \
        spectrum_file=xs.dat \
        spectun=0 \
        trad=0 \
        rlrad38=3.6545e-013 \
        column=3.0599e+016 \
        rlogxi=3.1217e+000 \
        lcpres=0 \
        Habund=1.000e+000 \
        Heabund=1.000e+000 \
        Liabund=0.000e+000 \
        Beabund=0.000e+000 \
        Babund=0.000e+000 \
        Cabund=0.000e+000 \
        Nabund=0.000e+000 \
        Oabund=0.000e+000 \
        Fabund=0.000e+000 \
        Neabund=0.000e+000 \
        Naabund=0.000e+000 \
        Mgabund=0.000e+000 \
        Alabund=0.000e+000 \
        Siabund=2.857e+005 \
        Pabund=0.000e+000 \
        Sabund=0.000e+000 \
        Clabund=0.000e+000 \
        Arabund=0.000e+000 \
        Kabund=0.000e+000 \
        Caabund=0.000e+000 \
        Scabund=0.000e+000 \
        Tiabund=0.000e+000 \
        Vabund=0.000e+000 \
        Crabund=0.000e+000 \
        Mnabund=0.000e+000 \
        Feabund=0.000e+000 \
        Coabund=0.000e+000 \
        Niabund=0.000e+000 \
        Cuabund=0.000e+000 \
        Znabund=0.000e+000 \
        modelname=run3 \
        niter=-99 \
        npass=3 \
        critf=1E-008 \
        nsteps=1 \
        xeemin=0.1 \
        emult=0.1 \
        taumax=20 \
        lwrite=1 \
        lprint=1 \
        lstep=0 \
        ncn2=999 \
        radexp=0 \
        vturbi=0')
        
os.system('chmod +rxw run_xstar.sh')
os.system('./run_xstar.sh')
os.system('mkdir '+directory_name)

output_files = ['xo01_detail.fits '
'xo01_detal2.fits '
'xo01_detal3.fits '
'xo01_detal4.fits '
'xo02_detail.fits '
'xo02_detal2.fits '
'xo02_detal3.fits '
'xo02_detal4.fits '
'xo03_detail.fits '
'xo03_detal2.fits '
'xo03_detal3.fits '
'xo03_detal4.fits '
'xout_abund1.fits '
'xout_cont1.fits '
'xout_lines1.fits '
'xout_rrc1.fits '
'xout_spect1.fits '
'xout_step.log']

for file in output_files:
    os.system('mv '+file+' '+directory_name)