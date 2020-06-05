#!/bin/csh -f
# run charmm test cases

#This script can take two optional variables, the first one to specify which charmm to use (default is the serial charmm on clusters) and the second one to specify a particulr testcase (default is to run all of them). After succesfully running the test, one should compare the resulted test_drude_all.ene with test_drude_all.ene.benchmark.

#charmm executable
#setenv charmm /home/alex/charmm/exec/charmm
#setenv charmm /home/huangj/charmm/charmm.drudefix/charmm

if ($#argv == 0) then
    echo "No charmm is specified, using default serial charmm on the cluster"
#    setenv charmm /opt/mackerell/apps/charmm/serial/charmm-serial
#    setenv charmm /opt/mackerell/apps/charmm/serial/c44b2-serial
    set all = 0
else if ($#argv == 1) then
    # need to make the judgement whether the user wants to specify charmm or specify a particular test to do
    if ( `echo $1 | grep '^[0-9]*$'`) then
#        setenv charmm /opt/mackerell/apps/charmm/serial/charmm-serial
#        setenv charmm /opt/mackerell/apps/charmm/serial/c44b2-serial
         set all = $1
    else
        # it's not a number, so should be charmm position
        setenv charmm $1
        set all = 0
    endif
else
    setenv charmm $1
    set all = $2
endif

setenv outdir  .

#mkdir -p $outdir

$charmm mindr:0 mini:0 test:$all < test_drude_all_2019_auto_all.inp > $outdir/test_drude_all_2019_auto_all_c44b2.out 

#grep for crashed jobs
grep '! I I I I I !' $outdir/*.out
grep '\   XXX   /' $outdir/*.out
grep 'CHECK THOSE INPUTS, ACE' $outdir/*.out
grep 'ABNORMAL TERMINATION' $outdir/*.out
