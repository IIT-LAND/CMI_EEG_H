# _03_cmi_postprocessing_hurst.py
"""
_03_cmi_postprocessing_hurst.py

Calculate the Hurst exponent on preprocessed CMI EEG data.

Example usage:

# Run on server and generate scripts only
python _03_cmi_postprocessing_hurst.py server script

# Run on laptop and generate scripts only
python _03_cmi_postprocessing_hurst.py laptop script

# Run on server and run batches
python _03_cmi_postprocessing_hurst.py server batch

# Run on laptop and run batches
python _03_cmi_postprocessing_hurst.py laptop batch

--- Written by nbertelsen and mvlombardo
"""

# import modules ---------------------------------------------------------------
import os
import sys
import fnmatch
import glob

# Parse input arguments --------------------------------------------------------
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

# specify number of computational threads for use n MATLAB
maxNumCompThreads = 1

metric2use = "H"
main_script = "cmi_postproc_hurst"

# Directories ------------------------------------------------------------------
# specify my rootdir that all other directories follow from
if run_on_server:
    rootdir = "/media/DATA/RAW/cmihbn"
else:
    rootdir = "/Users/mlombardo/Dropbox/cmi_test"

# other directories
codedir = "%s/code" % rootdir
preprocdir = "%s/data/preproc" % rootdir
postprocdir = "%s/data/postproc" % rootdir
os.system("mkdir -p %s" % postprocdir)


# Function to check if postproc was already done -------------------------------
def check_postproc_done(path2use, metric2use):
    # Function to check inside postproc directory if *metric2use*.csv file exists
    # to indicate that preprocessing was already run.

    fname = glob.glob("%s/*_%s.csv" % (path2use, metric2use))

    if not fname:
        postproc_complete = False
    else:
        postproc_complete = True

    return(postproc_complete)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

fstem1 = "_preproc_icaDenoised"
fstem2 = "_preproc_notDenoised"

# generate my list of subjects from subjects already present in the raw data directory
sublist = os.listdir(preprocdir)
# if not run_on_server:
#     sublist.remove('.DS_Store')

for sub in sublist:

    # subject-specific paths in preproc directories
    preproc_subpath = "%s/%s" % (preprocdir, sub)
    postproc_subpath = "%s/%s" % (postprocdir, sub)

    # generate a list of files inside subpath
    flist = os.listdir(preproc_subpath)
    # if not run_on_server:
    #     flist.remove('.DS_Store')
    # find other *.sh scripts in flist if the exist, and remove them
    pattern2find = "*.sh"
    matches = fnmatch.filter(flist, pattern2find)
    for match in matches:
        flist.remove(match)

    # loop over tasks in flist
    for task in flist:

        # check in postproc dir if *_metric2use*.set file exists to indicate
        # postprocessing was already done
        path2use = "%s/%s/%s" % (postproc_subpath,task,metric2use)
        postproc_complete = check_postproc_done(path2use, metric2use)

        if postproc_complete:
            # if postprocessing was already done, print this and move on...
            print("%s postprocessing already complete" % path2use)

        elif not postproc_complete:
            # else if postproc was not complete, run all the following...

            # start my bash script - create list
            script = []
            script.append("#!/bin/bash")
            script.append("")
            script.append("start=`date +%s`")
            script.append("")
            script.append("# %s %s" % (sub,metric2use))
            script.append("")
            script.append("# %s" % task)

            # make a task directory in the preprocessing directory
            indir = "%s/%s" % (preproc_subpath, task)
            outpath = "%s/%s/%s" % (postproc_subpath, task, metric2use)
            os.system("mkdir -p %s" % outpath)

            # file to compute the Hurst on
            datafile1 = "%s/%s_%s%s.set" % (indir,sub,task,fstem1)
            datafile2 = "%s/%s_%s%s.set" % (indir,sub,task,fstem2)
            if os.path.isfile(datafile1):
                datafile = datafile1
            elif os.path.isfile(datafile2):
                datafile = datafile2
            else:
                datafile = datafile1

            outfile = "%s/%s_%s_%s.csv" % (outpath, sub, task, metric2use)

            # generate code for the MATLAB function to run
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; [results, outfile] = %s('%s', true); cmi_topoplot(outfile); catch ME; disp(ME.message); toc; exit; end; toc; exit;" % (codedir, maxNumCompThreads, main_script, datafile)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)
            script.append("")
            script.append("end=`date +%s`")
            script.append("runtime=$((end-start))")
            script.append("echo Time to complete = $runtime secs")
            script.append("")

            # Writing bash script to file ------------------------------------------
            print("Saving postprocessing %s bash script for %s %s ..." %(metric2use,sub, task))
            fname2save = "%s/_03_postproc_%s_%s_%s.sh" % (outpath, metric2use,sub, task)
            fname = open(fname2save,"w")
            fname.write("\n".join(script)+"\n")
            fname.close()

            os.system("mkdir -p %s/_03_postproc_%s_batches" % (codedir,metric2use))
            fname2save = "%s/_03_postproc_%s_batches/_03_postproc_%s_%s_%s.sh" % (codedir,metric2use,metric2use,sub, task)
            fname = open(fname2save,"w")
            fname.write("\n".join(script)+"\n")
            fname.close()

            # if you're not just generating scripts and want to run the bash scripts
            # then do this below...
            if not script_only:

                if run_on_server:

                    print("Inserting batch into the queue ...")

                    # options for sbatch
                    outlog = "%s/%s_%s_postproc_%s_out.log" % (outpath, sub, task, metric2use)
                    errlog = "%s/%s_%s_postproc_%s_error.log" % (outpath, sub, task, metric2use)
                    sbatch_opts = "-J %s --tasks-per-node=1 --cpus-per-task=1 --time=06:00:00 --no-requeue --mem=8g --output=%s --error=%s" % (sub, outlog, errlog)

                    # run batch
                    os.system("sbatch %s %s" % (sbatch_opts, fname2save))

                else:

                    # run bash script
                    os.system("bash %s" % fname2save)
