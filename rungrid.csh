#!/bin/csh

umask 002 #Set permissions

CPN=/grid/fermiapp/minos/scripts/cpn
MVN=/grid/fermiapp/minos/scripts/mvn

mkdir ${_CONDOR_SCRATCH_DIR}/input 
mkdir ${_CONDOR_SCRATCH_DIR}/output

${CPN} /minos/data/analysis/AtmosNu/atnu_ntuples_dogwood/ntpsummary.mc.upmu.0000.root/ ${_CONDOR_SCRATCH_DIR}/input #Copy simulations to directory on grid machine

${CPN} /minos/data/analysis/AtmosNu/atnu_ntuples_dogwood/ntpsummary.mc.atmos.0000.root/ ${_CONDOR_SCRATCH_DIR}/input

# Copy data to directory on grid machines
${CPN} /minos/app/tjyang/atnu/barrflux/fmax20_0401z.sou_nue ${_CONDOR_SCRATCH_DIR}/input
${CPN} /minos/app/tjyang/atnu/barrflux/fmax20_0401z.sou_nbe ${_CONDOR_SCRATCH_DIR}/input
${CPN} /minos/app/tjyang/atnu/barrflux/fmax20_0401z.sou_num ${_CONDOR_SCRATCH_DIR}/input
${CPN} /minos/app/tjyang/atnu/barrflux/fmax20_0401z.sou_nbm ${_CONDOR_SCRATCH_DIR}/input

${CPN} /afs/fnal.gov/files/home/room3/corwin/MINOS/analyses/atmosnu/ana.* ${_CONDOR_SCRATCH_DIR}/input

export INPUTFILES=${_CONDOR_SCRATCH_DIR}/input/
export OUTPUTFILE=${_CONDOR_SCRATCH_DIR}/output/

loon -qb rungrid.C(0,1.0,2.43E-3,"$INPUTFILES","$OUTPUTFILE") #Run the analysis code

# Move files from grid machine back to Fermilab cluster
${MVN} $OUTPUTFILE/*.root /minos/data/users/corwin/masshierachy/ 
