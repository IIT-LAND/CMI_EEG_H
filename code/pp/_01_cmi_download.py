# _01_cmi_download.py
"""
_01_cmi_download.py

Writes and runs bash scripts for downloading subjects from the CMI dataset.

First argument should be either server or laptop to indicate to the script if you're running the script on the lab's server or your own laptop (for testing)

Second argument is which group to run it on. Should be either "TD", "ASD", "ADHD", "ID", or "LANG".

Third argument is either "script" or "batch", and will indicate if you want to generate scripts only or if you want to generate the scripts and then run the batches.


Example usage:

# Run on the server for TD for R1 to R10 and generate scripts only
python _01_cmi_download.py server TD script 1 10

# Run on the server for ASD for R1 to R10 and generate scripts AND run batches all in one go
python _01_cmi_download.py laptop ASD batch 1 10

--- Written by mvlombardo
"""

# import modules ---------------------------------------------------------------
import os
import sys
import numpy as np
import pandas as pd


# Parse input arguments --------------------------------------------------------
# First argument tells us whether we are running it on the lab's server or our own laptop
if sys.argv[1]=="server":
    run_on_server = True
elif sys.argv[1]=="laptop":
    run_on_server = False

# Second argument tells us which group to run on
if sys.argv[2]=="TD":
    group2run = "TD"
elif sys.argv[2]=="ASD":
    group2run = "ASD"
elif sys.argv[2]=="ADHD":
    group2run = "ADHD"
elif sys.argv[2]=="ID":
    group2run = "ID"
elif sys.argv[2]=="LANG":
    group2run = "LANG"

# Third argument tells if you are running to generate scripts only or to run batches as well
if sys.argv[3]=="script":
    script_only = True
elif sys.argv[3]=="batch":
    script_only = False

# Fourth argument tells the starting Data Release to use
start_data_release = np.array(sys.argv[4], dtype = int)

# Fifth argument tells the ending Data Release to use
end_data_release = np.array(sys.argv[5], dtype = int)

# specify number of computational threads for use n MATLAB
maxNumCompThreads = 1

# Directories ------------------------------------------------------------------
if run_on_server:
    rootdir="/media/DATA/RAW/cmihbn"
else:
    rootdir="/Users/mlombardo/Dropbox/cmi_test"


datadir="%s/data/raw" % rootdir
codedir="%s/code" % rootdir
roottmpdir="%s/tmp" % rootdir
tmpdir="%s/tmp" % datadir

os.system("mkdir -p %s" % roottmpdir)
os.system("mkdir -p %s" % tmpdir)


# Generate subject lists to loop over ------------------------------------------
os.system("Rscript %s/cmi_wrangle_pheno_data.R %s %d %d" % (codedir,group2run,start_data_release,end_data_release))

if start_data_release==end_data_release:
    fstem = "R%d" % start_data_release
else:
    fstem = "R%d_to_R%d" % (start_data_release,end_data_release)

if group2run=="TD":
    sublist2use = "%s/td_sublist_%s.csv" % (roottmpdir,fstem)
elif group2run=="ASD":
    sublist2use = "%s/asd_sublist_%s.csv" % (roottmpdir,fstem)
elif group2run=="ADHD":
    sublist2use = "%s/adhd_sublist_%s.csv" % (roottmpdir,fstem)
elif group2run=="ID":
    sublist2use = "%s/id_sublist_%s.csv" % (roottmpdir,fstem)
elif group2run=="LANG":
    sublist2use = "%s/lang_sublist_%s.csv" % (roottmpdir,fstem)

sublist = pd.read_csv(sublist2use, header = None)
sublist = sublist[0].tolist()

# Loop over subjects and write out their script --------------------------------
for subid in sublist:

    print(subid)

    # make a bash script for a specific subject --------------------------------
    script = []
    script.append("#!/bin/bash")

    # set up directories in the bash script ------------------------------------
    script.append("")
    script.append("# Setup directories")
    script.append("rootdir=%s" % rootdir)
    script.append("datadir=%s" % datadir)
    script.append("codedir=%s" % codedir)
    script.append("roottmpdir=%s" % roottmpdir)
    script.append("tmpdir=%s" % tmpdir)
    script.append("cd %s" % codedir)

    # set up subject-specific directory ----------------------------------------
    script.append("")
    script.append("# Subject info")
    script.append("subid=%s" % subid)
    subdir = "%s/%s" % (datadir,subid)
    script.append("subdir=%s" % (subdir))

    # download the data from the Amazon cloud ----------------------------------
    script.append("")
    script.append("# Download the data")
    cloud_url = "https://fcp-indi.s3.amazonaws.com/data/Archives/HBN/EEG/%s.tar.gz" % subid
    downloaded_tar_file = "%s/%s.tar.gz" % (tmpdir,subid)
    if run_on_server:
        script.append("cd %s" % tmpdir)
        script.append("wget %s" % cloud_url)
        script.append("cd %s" % codedir)
    elif not run_on_server:
        script.append("curl %s -o %s" % (cloud_url,downloaded_tar_file))

    # unpack the downloaded tar file -------------------------------------------
    script.append("")
    script.append("# Unpack the downloaded tar ball")
    script.append("cd %s" % tmpdir)
    script.append("tar -xzf %s" % downloaded_tar_file)
    script.append("cd %s" % codedir)
    script.append("mkdir -p %s" % subdir)

    # RestingState -------------------------------------------------------------
    dataType = "RestingState"
    script.append("")
    script.append("# %s" % dataType)
    # check if .mat file exists
    fname2find = "%s/%s/EEG/raw/mat_format/%s.mat" % (tmpdir,subid,dataType)
    script.append('if test -f "%s"; then' % fname2find)
    script.append('    destdir=%s/%s' % (subdir,dataType))
    dest_dir = "%s/%s" % (subdir,dataType)
    script.append('    mkdir -p ${destdir}')
    newfname = "%s/%s_%s.mat" % (dest_dir,subid,dataType)
    script.append('    cp -Rf %s %s' % (fname2find,newfname))
    script.append("fi")
    # # --------------------------------------------------------------------------

    # Video1 -------------------------------------------------------------------
    dataType = "Video1"
    script.append("")
    script.append("# %s" % dataType)
    # check if .mat file exists
    fname2find = "%s/%s/EEG/raw/mat_format/%s.mat" % (tmpdir,subid,dataType)
    script.append('if test -f "%s"; then' % fname2find)
    script.append('    destdir=%s/%s' % (subdir,dataType))
    dest_dir = "%s/%s" % (subdir,dataType)
    script.append('    mkdir -p ${destdir}')
    newfname = "%s/%s_%s.mat" % (dest_dir,subid,dataType)
    script.append('    cp -Rf %s %s' % (fname2find,newfname))
    script.append('    newfname=%s' % newfname)
    matlab_commands2run = "cd('${codedir}'); maxNumCompThreads(%d); cmi_rename_videos('${newfname}'); exit;" % maxNumCompThreads
    script.append('    matlab -nodesktop -r "%s"' % matlab_commands2run)
    script.append("fi")

    # Video2 -------------------------------------------------------------------
    dataType = "Video2"
    script.append("")
    script.append("# %s" % dataType)
    # check if .mat file exists
    fname2find = "%s/%s/EEG/raw/mat_format/%s.mat" % (tmpdir,subid,dataType)
    script.append('if test -f "%s"; then' % fname2find)
    script.append('    destdir=%s/%s' % (subdir,dataType))
    dest_dir = "%s/%s" % (subdir,dataType)
    script.append('    mkdir -p ${destdir}')
    newfname = "%s/%s_%s.mat" % (dest_dir,subid,dataType)
    script.append('    cp -Rf %s %s' % (fname2find,newfname))
    script.append('    newfname=%s' % newfname)
    matlab_commands2run = "cd('${codedir}'); maxNumCompThreads(%d); cmi_rename_videos('${newfname}'); exit;" % maxNumCompThreads
    script.append('    matlab -nodesktop -r "%s"' % matlab_commands2run)
    script.append("fi")

    # Video3 -------------------------------------------------------------------
    dataType = "Video3"
    script.append("")
    script.append("# %s" % dataType)
    # check if .mat file exists
    fname2find = "%s/%s/EEG/raw/mat_format/%s.mat" % (tmpdir,subid,dataType)
    script.append('if test -f "%s"; then' % fname2find)
    script.append('    destdir=%s/%s' % (subdir,dataType))
    dest_dir = "%s/%s" % (subdir,dataType)
    script.append('    mkdir -p ${destdir}')
    newfname = "%s/%s_%s.mat" % (dest_dir,subid,dataType)
    script.append('    cp -Rf %s %s' % (fname2find,newfname))
    script.append('    newfname=%s' % newfname)
    matlab_commands2run = "cd('${codedir}'); maxNumCompThreads(%d); cmi_rename_videos('${newfname}'); exit;" % maxNumCompThreads
    script.append('    matlab -nodesktop -r "%s"' % matlab_commands2run)
    script.append("fi")

    # Video4 -------------------------------------------------------------------
    dataType = "Video4"
    script.append("")
    script.append("# %s" % dataType)
    # check if .mat file exists
    fname2find = "%s/%s/EEG/raw/mat_format/%s.mat" % (tmpdir,subid,dataType)
    script.append('if test -f "%s"; then' % fname2find)
    script.append('    destdir=%s/%s' % (subdir,dataType))
    dest_dir = "%s/%s" % (subdir,dataType)
    script.append('    mkdir -p ${destdir}')
    newfname = "%s/%s_%s.mat" % (dest_dir,subid,dataType)
    script.append('    cp -Rf %s %s' % (fname2find,newfname))
    script.append('    newfname=%s' % newfname)
    matlab_commands2run = "cd('${codedir}'); maxNumCompThreads(%d); cmi_rename_videos('${newfname}'); exit;" % maxNumCompThreads
    script.append('    matlab -nodesktop -r "%s"' % matlab_commands2run)
    script.append("fi")

    # Remove the temporary downloaded data -------------------------------------
    script.append("")
    script.append("# Remove temporary downloaded data")
    script.append("rm -Rf %s/%s*" % (tmpdir,subid))
    script.append("# Done!")

    # Writing bash script to file ----------------------------------------------
    print("Saving bash script ...")

    fname2save = "%s/_01_download_batches/_01_download_%s.sh" % (codedir,subid)
    os.system("mkdir -p %s/_01_download_batches" % codedir)
    fname = open(fname2save,"w")
    fname.write("\n".join(script)+"\n")
    fname.close()

    # Run batches --------------------------------------------------------------
    if not script_only:

        if run_on_server:
            print("Inserting batch into the queue ...")

            # options for sbatch
            sbatch_opts = '-J %s --tasks-per-node=1 --cpus-per-task=1 --time=06:00:00 --no-requeue --mem=8g' % subid

            # run batch
            os.system("sbatch %s %s" % (sbatch_opts,fname2save))

        else:

            # run bash script
            os.system("bash %s" % fname2save)

print("Done!")
