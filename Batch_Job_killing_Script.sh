bjobs | awk {'print "bkill ",$1'} > Kill_Batch_Jobs.sh
source Kill_Batch_Jobs.sh
rm Kill_Batch_Jobs.sh
