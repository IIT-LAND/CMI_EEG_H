# _04_cmi_tidy_hurst.py
"""
_04_cmi_tidy_hurst.py

Put all subjects Hurst exponent data into a tidy csv file.

Example usage:

# Run on server and generate scripts only
python _04_cmi_tidy_hurst.py server script

# Run on laptop and generate scripts only
python _04_cmi_tidy_hurst.py laptop script

# Run on server and run batches
python _04_cmi_tidy_hurst.py server batch

# Run on laptop and run batches
python _04_cmi_tidy_hurst.py laptop batch

--- Written by mvlombardo
"""

# import modules
import os
import sys
import fnmatch
import glob
import time

# Parse input arguments
# First argument tells us whether we are running it on the lab's server or our own laptop
if sys.argv[1]=="server":
    run_on_server = True
elif sys.argv[1]=="laptop":
    run_on_server = False

# Second argument tells if you are running to generate scripts only or to run batches as well
if sys.argv[2]=="script":
    script_only = True
elif sys.argv[2]=="batch":
    script_only = False

metric2use = "H"
main_script = "cmi_tidy_hurst"
# timestr = time.strftime("%Y%m%d-%H%M%S")
timestr = time.strftime("%Y%m%d")

# Directories
# specify my rootdir that all other directories follow from
if run_on_server:
    rootdir = "/media/DATA/RAW/cmihbn"
else:
    rootdir = "/Users/mlombardo/Dropbox/cmi_test"

# other directories
codedir = "%s/code" % rootdir
tmpdir = "%s/data/tmp" % rootdir
preprocdir = "%s/data/preproc" % rootdir
postprocdir = "%s/data/postproc" % rootdir
tidydir = "%s/data/tidy" % rootdir
os.system("mkdir -p %s" % tidydir)

# Initialize masterfiles and sublists
masterfiles = {}
sublists = {}
task_names = ["RestingState","Despicable","Fractals","Wimpy","Present"]
for task in task_names:
    # masterfile names
    masterfiles[task] = "%s/tidy_%s_%s_%s.csv" % (tidydir,task,metric2use,timestr)
    # initialize sublists
    sublists[task] = []

# generate my list of subjects from subjects already present in the raw data directory
sublist = os.listdir(postprocdir)
# if not run_on_server:
#     sublist.remove('.DS_Store')

# loop over subjects
for sub in sublist:

    # subject-specific paths in preproc directories
    postproc_subpath = "%s/%s" % (postprocdir, sub)

    # generate a list of files inside subpath
    flist = os.listdir(postproc_subpath)
    # find other *.sh scripts in flist if the exist, and remove them
    pattern2find = "*.sh"
    matches = fnmatch.filter(flist, pattern2find)
    for match in matches:
        flist.remove(match)

    # loop over tasks in flist
    for task in flist:
        # append subject to task-specific subject list

        taskpath = "%s/%s/%s" % (postproc_subpath, task, metric2use)
        file2lookfor = "%s/%s_%s_%s.csv" %(taskpath,sub,task,metric2use)
        if os.path.isfile(file2lookfor):
            sublists[task].append(sub)


# Writing out batch file
outpath = "%s/_04_tidy_%s_batches" % (codedir,metric2use)
os.system("mkdir -p %s" % outpath)

# start my bash script - create list
script = []
script.append("#!/bin/bash")
script.append("")
script.append("start=`date +%s`")
script.append("")

sublist2save = {}
for task in task_names:

    # save out temporary sublist
    sublist2save[task] = "%s/_04_tidy_%s_%s_%s.txt" % (outpath, metric2use, task, timestr)
    fname = open(sublist2save[task],"w")
    fname.write("\n".join(sublists[task])+"\n")
    fname.close()

    # R markdown to make the report
    script.append("")
    r_code2run = "source('%s/cmi_tidy_hurst.R'); task = '%s'; sublist = read.delim('%s',header=FALSE); masterfile = '%s'; result_df = cmi_tidy_hurst(sublist,task, masterfile); quit(save = 'no', status = 0, runLast = FALSE);" % (codedir,task, sublist2save[task],masterfiles[task])
    script.append('Rscript -e "%s"' % r_code2run)

script.append("")
script.append("end=`date +%s`")
script.append("runtime=$((end-start))")
script.append("echo Time to complete = $runtime secs")
script.append("")

# write batch out to file
fname2save = "%s/_04_tidy_%s_%s.sh" % (outpath,metric2use,timestr)
fname = open(fname2save,"w")
fname.write("\n".join(script)+"\n")
fname.close()

# if you're not just generating scripts and want to run the bash scripts
# then do this below...
if not script_only:

    if run_on_server:

        print("Inserting batch into the queue ...")

        # options for sbatch
        outlog = "%s/_04_tidy_%s_%s_out.log" % (outpath, metric2use,timestr)
        errlog = "%s/_04_tidy_%s_%s_error.log" % (outpath, metric2use,timestr)
        sbatch_opts = "-J tidy_%s --tasks-per-node=1 --cpus-per-task=1 --time=99:00:00 --no-requeue --mem=8g --output=%s --error=%s" % (metric2use, outlog, errlog)

        # run batch
        os.system("sbatch %s %s" % (sbatch_opts, fname2save))

    else:

        # run bash script
        os.system("bash %s" % fname2save)
