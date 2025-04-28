#!/usr/bin/env python
"""
-----------------------------------------------------------------------
 Description:  Example script to submit a job through the HPC batch system

 Assumptions: For use on Hera. This script must be copied to the bin directory.
              The user should edit the JOB_NAME, ACCOUNT, etc.
              ./scm/etc/Hera_setup_intel.sh or
              ./scm/etc/Hera_setup_intel.csh must be run before submitting this script

 COMMAND line arguments: none

 Usage: ./scm_slurm_example.py   # from the bin directory
-----------------------------------------------------------------------
"""

from subprocess import Popen, PIPE
import os
import time

USER = os.getenv('USER')
STRINGS = [USER, 'noaa.gov']
MY_EMAIL = '@'.join(STRINGS)
PROC = Popen('sbatch', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

### User-editable section ###
JOB_NAME = "scmout"
ACCOUNT = "ufs-phys"
WALLTIME = "20"
PROCESSORS = "1"
QUEUE = "batch"
#COMMAND = "./run_scm.py --npz_type gfs --file scm_ufsens_control_c192.py --timestep 360"
#COMMAND = "./run_scm.py --npz_type gfs --file scm_ufsens_HYD_control_p8_ugwpv1_intel.py --timestep 600 --n_itt_out 6 --n_itt_diag 6"
#COMMAND = "./run_scm.py --npz_type gfs --file scm_ufsens_control_p8_ugwpv1_intel.py --timestep 600 --n_itt_out 6 --n_itt_diag 6"
COMMAND = "./cmp_scmout.py -fbl /scratch2/BMC/ufs-phys/Man.Zhang/_case_gen/ccpp-scm_v7_hr_v4/scm/run/HR4_C96_NH/output_control_p8_ugwpv1_intel_n001_SCM_GFS_v17_p8_ugwpv1/output.nc -frt /scratch2/BMC/ufs-phys/Man.Zhang/_case_gen/ccpp-scm_v7_hr_v4/scm/run/HR4_C96_HYD/output_HYD_control_p8_ugwpv1_intel_n001_SCM_GFS_v17_p8_ugwpv1/output.nc"
EMAIL_ADDR = MY_EMAIL
EMAIL_OCCASION = "BEGIN,END,FAIL"
SERIAL_MEM = "512M"
WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
### End User-editable section ###


JOB_STRING = """#!/bin/bash
#SBATCH -J %s
#SBATCH -A %s
#SBATCH --time=%s
#SBATCH -n %s
#SBATCH -o ./%s.out
#SBATCH -q %s
#SBATCH --mail-user=%s
#SBATCH --mail-type=%s
#SBATCH --mem=%s
#SBATCH -D%s
%s""" % (JOB_NAME, ACCOUNT, WALLTIME, PROCESSORS, JOB_NAME, QUEUE, EMAIL_ADDR,
         EMAIL_OCCASION, SERIAL_MEM, WORKING_DIR, COMMAND)


# Send JOB_STRING to qsub
PROC.stdin.write(JOB_STRING.encode())
OUT, ERR = PROC.communicate()

# Print your job and the system response to the screen as it's submitted
print(JOB_STRING)
print(OUT.decode())
print(ERR.decode())

time.sleep(0.1)
