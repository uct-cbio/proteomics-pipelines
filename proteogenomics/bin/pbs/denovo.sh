#PBS -N MyJob
#PBS -q UCTlong
#PBS -l nodes=1:ppn=20:series600
#PBS -M matthys.potgieter@gmail.com
#PBS -m ae

# NB, PLEASE READ THIS!
# There is a 2:1 correspondence between RAM and cores on the 600 series.
# You need to know how much RAM your job will consume before submitting it.
# Please set the ppn value above to be 1/2 the GB of RAM required.  For
# example a job needing 10GB of RAM should have ppn=5

# Please leave the hostname command here for troubleshooting purposes.
hostname

# Your science stuff goes here:
DeNovoGUI.sh /home/ptgmat003/cbio-pipelines/proteogenomics/bin/config/DeNovoGui_example.cfg



