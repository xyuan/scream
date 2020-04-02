#!/bin/sh

testname=`date +"%Y_%m_%d_%H%M%S"`
echo "Running test for date: $testname"
./create_test \
  SMS_P36x1.ne4_ne4.FSCREAM-SA \
  --machine quartz \
  --baseline-root  /p/lustre2/donahue5/E3SM_simulations/SCREAM/scream_development/SCORPIO_ingegrations/\
  -r /p/lustre2/donahue5/E3SM_simulations/SCREAM/scream_development/SCORPIO_ingegrations/ \
  --output-root /p/lustre2/donahue5/E3SM_simulations/SCREAM/scream_development/SCORPIO_ingegrations/ \
  --project cbronze \
  -q pdebug \
  -t $testname

#  ERS_P36x2.ne4_ne4.FC5AV1C-L \
#  ERS_P48x1.ne4_ne4.FC5AV1C-L \
#  ERS_P96x1.ne4_ne4.FC5AV1C-L \
#  ERS_P24x2.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
#  ERS_P48x2.ne4_ne4.FSCREAM-LR.quartz_intel.cam-force_smp \
