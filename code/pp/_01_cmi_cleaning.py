# _01_cmi_cleaning.py
"""
_01_cmi_cleaning.py

Goes through the initial downloaded data and does the following:

1) Step 1:
        Moves empty directories to a temporary directory (because data wasn't
            downloaded for whatever reason).
        Moves directories of videos that couldn't be renamed to a temporary
            directory.
2) Step 2:
        Insert channel info into EEG structure.
        Removes irrelevant data segments.
        Adds latencies to EEG structure.
        Saves to a *.set file.

3) Step 3:
        Removes all Video* directories that are still left over...

Steps 1-3 are run separately, which means calling the script each time
changing the second argument to either step1, step2, or step3

Example usage:

# Run on laptop, step1 and generate scripts only
python _01_cmi_cleaning.py laptop step1 script

# Run on the server, step2, AND run batches all in one go
python _01_cmi_cleaning.py server step2 batch

# Run all steps on the server in batch modules
python _01_cmi_cleaning.py server step1 batch
python _01_cmi_cleaning.py server step2 batch
python _01_cmi_cleaning.py server step3 batch

--- Written by mvlombardo
"""

# import modules ---------------------------------------------------------------
import os
import sys
import fnmatch

# ------------------------------------------------------------------------------
# function to check if a Video directory exists
def check_for_video_dir(flist):

    video_names = list()
    for fname in flist:
        if "Video" in fname:
            video_names.append(fname)

    if any("Video" in fname for fname in flist):
        video_dir_exists = True
    else:
        video_dir_exists = False

    return(video_dir_exists,video_names)

# Parse input arguments --------------------------------------------------------
# First argument tells us whether we are running it on the lab's server or our own laptop
if sys.argv[1]=="server":
    run_on_server = True
elif sys.argv[1]=="laptop":
    run_on_server = False

# Second argument tells which step you are running
if sys.argv[2]=="step1":
    run_step1 = True
    run_step2 = False
    run_step3 = False
elif sys.argv[2]=="step2":
    run_step1 = False
    run_step2 = True
    run_step3 = False
elif sys.argv[2]=="step3":
    run_step1 = False
    run_step2 = False
    run_step3 = True

# Third argument tells if you are running to generate scripts only or to run batches as well
if sys.argv[3]=="script":
    script_only = True
elif sys.argv[3]=="batch":
    script_only = False

# specify number of computational threads for use n MATLAB
maxNumCompThreads = 1

# Directories ------------------------------------------------------------------
if run_on_server:
    rootdir = "/media/DATA/RAW/cmihbn"
else:
    rootdir = "/Users/mlombardo/Dropbox/cmi_test"

datadir = "%s/data/raw" % rootdir
codedir = "%s/code" % rootdir
roottmpdir = "%s/tmp" % rootdir
tmpdir = "%s/tmp" % datadir
problemdir = "%s/problematic_data" % tmpdir
os.system("mkdir -p %s" % problemdir)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Step 1:  Move all empty directories or Video data that wasn't renamed properly
# to a temporary directory
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if run_step1:

    # Subject list -------------------------------------------------------------
    sublist = os.listdir(datadir)
    # sublist.remove('.DS_Store')
    sublist.remove('tmp')

    script = []
    script.append("#!/bin/bash")

    for sub in sublist:

        subpath = "%s/%s" % (datadir,sub)
        flist = os.listdir(subpath)
        pattern2find = "*.sh"
        matches = fnmatch.filter(flist, pattern2find)
        for match in matches:
            flist.remove(match)

        (video_dir_exists, video_names) = check_for_video_dir(flist)

        if video_dir_exists:

            script.append("")
            script.append("mkdir -p %s/%s" %(problemdir,sub))
            script.append("mv %s/Video* %s/%s" % (subpath, problemdir,sub))
            error_message = "!!! Problem with data. Could not rename. !!!"

            for video_name in video_names:
                script.append("echo %s >> %s/%s/%s/renaming_error.txt" % (error_message,problemdir,sub,video_name))

        elif not flist:

            script.append("")
            script.append("mv %s %s" % (subpath, problemdir))
            error_message = "!!! No data downloaded !!!"
            script.append("echo %s >> %s/%s/download_error.txt" % (error_message,problemdir,sub))


    # Writing bash script to file ----------------------------------------------
    print("Saving cleaning bash script ...")
    fname2save = "%s/_01a_cleaning.sh" % (codedir)
    fname = open(fname2save,"w")
    fname.write("\n".join(script)+"\n")
    fname.close()
    if not script_only:
        os.system("bash %s" % fname2save)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Step 2:  Clean data to remove all irrelevant data segments and save as *.set
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if run_step2:

    # channel info file
    chan_info_file = "%s/GSN-HydroCel-129.sfp" % codedir

    # Subject list -------------------------------------------------------------
    sublist = os.listdir(datadir)
    # sublist.remove('.DS_Store')
    sublist.remove('tmp')

    # loop over subjects
    for sub in sublist:

        script = []
        script.append("#!/bin/bash")
        script.append("")
        script.append("# %s" % sub)

        # keep a running total of how many cleaning jobs are needed
        count_to_do = 0

        # make subpath
        subpath = "%s/%s" % (datadir,sub)
        # get a list of tasks
        flist = os.listdir(subpath)
        pattern2find = "*.sh"
        matches = fnmatch.filter(flist, pattern2find)
        for match in matches:
            flist.remove(match)

        # loop over tasks
        for task in flist:

            # make task path
            taskpath = "%s/%s" % (subpath,task)
            # make datafile
            datafile = "%s/%s_%s.mat" % (taskpath,sub,task)
            datafile2lookfor = "%s/%s_%s.set" % (taskpath,sub,task)

            if os.path.exists(datafile2lookfor):

                print("Skipping %s because file already exists" % datafile2lookfor)

            else:

                # update count_to_do
                count_to_do = count_to_do + 1

                # format call to MATLAB function to do cleaning
                script.append("")
                script.append("# Clean %s %s" % (sub,task))
                if run_on_server:
                    matlab_commands2run = "cd('%s'); maxNumCompThreads(%d); cmi_segment_clean('%s','%s',true); exit;" % (codedir,maxNumCompThreads,datafile,chan_info_file)
                else:
                    matlab_commands2run = "cd('%s'); maxNumCompThreads(%d); cmi_segment_clean('%s','%s',false); exit;" % (codedir,maxNumCompThreads,datafile,chan_info_file)

                script.append('matlab -nodesktop -r "%s"' % matlab_commands2run)

        # Writing bash script to file ------------------------------------------
        if count_to_do>0:
            print("Saving cleaning bash script ...")
            fname2save = "%s/_01b_cleaning_%s.sh" % (subpath,sub)
            fname = open(fname2save,"w")
            fname.write("\n".join(script)+"\n")
            fname.close()

            fname2save = "%s/_01b_cleaning_batches/_01b_cleaning_%s.sh" % (codedir,sub)
            os.system("mkdir -p %s/_01b_cleaning_batches" % codedir)
            fname = open(fname2save,"w")
            fname.write("\n".join(script)+"\n")
            fname.close()

            # Run batches ----------------------------------------------------------
            if not script_only:

                if run_on_server:
                    print("Inserting batch into the queue ...")

                    # options for sbatch
                    sbatch_opts = '-J %s --tasks-per-node=1 --cpus-per-task=1 --time=02:00:00 --no-requeue --mem=8g' % sub

                    # run batch
                    os.system("sbatch %s %s" % (sbatch_opts,fname2save))

                else:

                    # run bash script
                    os.system("bash %s" % fname2save)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Step 3:  Remove all Video* directories
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if run_step3:

    # Subject list -------------------------------------------------------------
    sublist = os.listdir(datadir)
    # sublist.remove('.DS_Store')
    sublist.remove('tmp')

    script = []
    script.append("#!/bin/bash")
    script.append("")

    for sub in sublist:

        subpath = "%s/%s" % (datadir,sub)
        flist = os.listdir(subpath)
        pattern2find = "*.sh"
        matches = fnmatch.filter(flist, pattern2find)
        for match in matches:
            flist.remove(match)

        for task in flist:

            path2check = "%s/%s" %(subpath, task)
            flist2 = os.listdir(path2check)
            (video_dir_exists, video_names) = check_for_video_dir(flist2)

            if video_dir_exists:
                script.append("")
                script.append("rm -Rf %s/Video*" % path2check)

    # Writing bash script to file ----------------------------------------------
    print("Saving cleaning bash script ...")
    fname2save = "%s/_01c_cleaning.sh" % (codedir)
    fname = open(fname2save,"w")
    fname.write("\n".join(script)+"\n")
    fname.close()
    if not script_only:
        os.system("bash %s" % fname2save)


print("Done!")
