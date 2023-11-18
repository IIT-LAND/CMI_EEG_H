# _02_cmi_preprocessing.py
"""
_02_cmi_preprocessing.py

Run preprocessing on CMI EEG data. 

Example usage:

# Run on server and generate scripts only
python _02_cmi_preprocessing.py server script step1

# Run on server and generate scripts and run in batch mode
python _02_cmi_preprocessing.py server batch step1

--- Written by mvlombardo
"""

# import modules ---------------------------------------------------------------
import os
import sys
import fnmatch
import glob
import time

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

# Third argument what step you are running
if sys.argv[3]=="step1":
    step1 = True
    check_preproc = False
    preproc_plots = False
elif sys.argv[3]=="check":
    step1 = False
    check_preproc = True
    preproc_plots = False
    timestr = time.strftime("%Y%m%d-%H%M%S")
elif sys.argv[3]=="plots":
    step1 = False
    check_preproc = False
    preproc_plots = True

# specify number of computational threads for use n MATLAB
maxNumCompThreads = 1

# Directories ------------------------------------------------------------------
# specify my rootdir that all other directories follow from
if run_on_server:
    rootdir = "/media/DATA/RAW/cmihbn"
    server_flag = "true"
else:
    rootdir = "/Users/mlombardo/Dropbox/cmi_test"
    server_flag = "false"

# other directories I need
datadir = "%s/data/raw" % rootdir
rawdir = "%s/data/raw" % rootdir
preprocdir = "%s/data/preproc" % rootdir
codedir = "%s/code" % rootdir
roottmpdir = "%s/tmp" % rootdir
tmpdir = "%s/tmp" % datadir
os.system("mkdir -p %s" % preprocdir)

# Function to check if preprocess was already done -----------------------------
def check_preproc_done(path2use):
    # Function to check inside preproc directory if *.nobadICA.set file exists
    # to indicate that preprocessing was already run.

    fname = glob.glob("%s/*Denoised.set" % path2use)

    if not fname:
        preproc_complete = False
    else:
        preproc_complete = True

    return(preproc_complete)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Main code
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

if step1:

    # generate my list of subjects from subjects already present in the raw data directory
    sublist = os.listdir(datadir)
    # sublist.remove('.DS_Store')
    sublist.remove('tmp')

    # preproc out filestem
    fstem1 = "_preproc_icaDenoised"
    fstem2 = "_preproc_notDenoised"

    # loop over subjects
    for sub in sublist:

        # subject-specific paths in raw and preproc directories
        subpath = "%s/%s" % (datadir,sub)
        pp_subpath = "%s/%s" % (preprocdir,sub)
        raw_subpath = "%s/%s" % (rawdir,sub)

        # generate a list of files inside subpath
        flist = os.listdir(subpath)

        # I want my flist to be just the tasks of interest. Other files like *.sh
        # files might exist in subpath, so here I need to find them and remove them
        # from flist before moving on. The remaining flist list is a list of strings
        # for each task directory
        pattern2find = "*.sh"
        matches = fnmatch.filter(flist, pattern2find)
        for match in matches:
            flist.remove(match)

        # loop over tasks in flist
        for task in flist:

            # check in preproc dir if *_nobadICA.set file exists to indicate
            # preprocessing was already done
            preproc_complete = check_preproc_done("%s/%s" % (pp_subpath,task))

            if preproc_complete:
                # if preprocessing was already done, print this and move on...
                print("%s %s preprocessing already complete" % (pp_subpath, task))

            elif not preproc_complete:
                # else if preproc was not complete, run all the following...

                # start my bash script
                script = []
                script.append("#!/bin/bash")
                script.append("")
                script.append("# %s Preprocessing" % sub)
                script.append("")
                script.append("# %s" % task)

                # make a task directory in the preprocessing directory
                os.system("mkdir -p %s/%s" % (pp_subpath,task))

                # file to preprocess
                datafile = "%s/%s/%s_%s.set" % (subpath,task,sub,task)

                # generate code for the MATLAB function to run
                function2run = "'cmi_run_preproc(\'\'%s\'\',%s)'" % (datafile,server_flag)

                # generate options to push through publish in MATLAB
                code2eval = "%s" % function2run
                options2use = "'format','pdf','outputDir','%s/%s', 'codeToEvaluate', %s" % (pp_subpath,task,code2eval)
                outlog = "%s/%s/%s_%s_preproc_matlab_output.txt" % (pp_subpath,task,sub,task)
                # generate the final MATLAB command to run
                matlab_command2run = "cd('%s'); diary '%s'; tic; maxNumCompThreads(%d); try; cmi_run_preproc('%s',%s); catch ME; disp(ME.message); toc; diary off; exit; end; toc; diary off; exit;" % (codedir,outlog,maxNumCompThreads,datafile,server_flag)
                script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

                script.append("")
                script.append("# %s Plotting and Report" % sub)
                script.append("")

                # file to preprocess
                outpath = "%s/%s" % (pp_subpath, task)
                outfile1 = "%s/%s_%s%s.set" % (outpath,sub,task,fstem1)
                outfile2 = "%s/%s_%s%s.set" % (outpath,sub,task,fstem2)
                rawpath2use = "%s/%s" % (raw_subpath,task)
                rawfile = "%s/%s_%s.set" % (rawpath2use, sub, task)

                # MATLAB stuff to for running scorepochs
                matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_scorepochs('%s'); catch; cmi_scorepochs('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,outfile1,outfile2)
                script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

                # MATLAB stuff to make the plots
                matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_tsplot('%s'); catch; cmi_tsplot('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,outfile1,outfile2)
                script.append('matlab -nodesktop -r "%s"' % matlab_command2run)
                matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_carpetplot('%s'); catch; cmi_carpetplot('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,outfile1,outfile2)
                script.append('matlab -nodesktop -r "%s"' % matlab_command2run)
                matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_tsnrplot('%s'); catch; cmi_tsnrplot('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,outfile1,outfile2)
                script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

                # MATLAB stuff to compute data quality metrics
                matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_preproc_quality('%s','%s'); catch; cmi_preproc_quality('%s','%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,rawfile,outfile1,rawfile,outfile2)
                script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

                # R markdown to make the report
                script.append("")
                r_code2run = "library(rmarkdown); setwd('%s'); render(input='%s/cmi_preproc_report.Rmd',output_file = '%s/%s_%s_preproc_report.html', params = list(subid = '%s', tasks = '%s'));" % (outpath, codedir,outpath,sub, task,sub, task)
                script.append('Rscript -e "%s"' % r_code2run)

                # Writing bash script to file ----------------------------------------------
                print("Saving preprocessing bash script for %s %s ..." %(sub,task))
                fname2save = "%s/%s/_02_preprocessing_%s_%s.sh" % (pp_subpath,task,sub,task)
                fname = open(fname2save,"w")
                fname.write("\n".join(script)+"\n")
                fname.close()

                os.system("mkdir -p %s/_02_preprocessing_batches" % codedir)
                fname2save = "%s/_02_preprocessing_batches/_02_preprocessing_%s_%s.sh" % (codedir,sub,task)
                fname = open(fname2save,"w")
                fname.write("\n".join(script)+"\n")
                fname.close()

                # if you're not just generating scripts and want to run the bash scripts
                # then do this below...
                if not script_only:

                    if run_on_server:

                        print("Inserting batch into the queue ...")

                        # options for sbatch
                        outlog = "%s/%s/%s_%s_preproc_out.log" % (pp_subpath,task,sub,task)
                        errlog = "%s/%s/%s_%s_preproc_error.log" % (pp_subpath,task,sub,task)
                        sbatch_opts = "-J %s --tasks-per-node=1 --cpus-per-task=1 --time=24:00:00 --no-requeue --mem=8g --output=%s --error=%s" % (sub,outlog,errlog)

                        # run batch
                        os.system("sbatch %s %s" % (sbatch_opts,fname2save))

                    else:

                        # run bash script
                        os.system("bash %s" % fname2save)


if check_preproc:

    # make the pp_issues_dir
    pp_issues_dir = "%s/data/tmp/preproc_issues_%s" % (rootdir,timestr)
    os.system("mkdir -p %s" % pp_issues_dir)
    os.system("mkdir -p %s/error_logs" % pp_issues_dir)
    os.system("mkdir -p %s/out_logs" % pp_issues_dir)

    # generate my list of subjects from preprocessed directory
    sublist = os.listdir(preprocdir)
    # sublist.remove('.DS_Store')
    # sublist.remove('tmp')

    # start my bash script
    script = []
    script.append("# Preprocessing check")
    script.append("")
    script.append("# Find the subjects and tasks where preprocessing did not complete because of issues")
    script.append("")
    script.append("# Check the pp_issues_dir for the *error.log and *out.log")
    script.append("# pp_issues_dir=%s" % pp_issues_dir)
    script.append("")

    issue_counter = 0
    # loop over subjects
    for sub in sublist:

        # subject-specific paths in raw and preproc directories
        pp_subpath = "%s/%s" % (preprocdir,sub)

        # generate a list of files inside subpath
        flist = os.listdir(pp_subpath)

        # I want my flist to be just the tasks of interest. Other files like *.sh
        # files might exist in subpath, so here I need to find them and remove them
        # from flist before moving on. The remaining flist list is a list of strings
        # for each task directory
        pattern2find = "*.sh"
        matches = fnmatch.filter(flist, pattern2find)
        for match in matches:
            flist.remove(match)

        # loop over tasks in flist
        for task in flist:

            # check in preproc dir if *_nobadICA.set file exists to indicate
            # preprocessing was already done
            preproc_complete = check_preproc_done("%s/%s" % (pp_subpath,task))

            if preproc_complete:
                # if preprocessing was already done, print this and move on...
                print("%s %s preprocessing already complete" % (pp_subpath, task))

            elif not preproc_complete:
                # else if preproc was not complete, run all the following...

                issue_counter = issue_counter+1
                path_of_interest = "%s/%s" % (pp_subpath, task)
                print("Preprocessing issue #%d found for %s" % (issue_counter, path_of_interest))
                script.append("%s" % path_of_interest)
                os.system("cp -Rf %s/*error.log %s/error_logs" % (path_of_interest,pp_issues_dir))
                os.system("cp -Rf %s/*out.log %s/out_logs" % (path_of_interest,pp_issues_dir))

    # Writing bash script to file ----------------------------------------------
    print("Saving pp_issues file ...")
    fname2save = "%s/_02_pp_issues_%s.txt" % (pp_issues_dir,timestr)
    fname = open(fname2save,"w")
    fname.write("\n".join(script)+"\n")
    fname.close()
    print("%s" % pp_issues_dir)


if preproc_plots:

    # preproc out filestem
    fstem1 = "_preproc_icaDenoised"
    fstem2 = "_preproc_notDenoised"

    # generate my list of subjects from subjects already present in the preproc data directory
    sublist = os.listdir(preprocdir)

    # loop over subjects
    for sub in sublist:

        # subject-specific paths in raw and preproc directories
        pp_subpath = "%s/%s" % (preprocdir,sub)
        raw_subpath = "%s/%s" % (rawdir,sub)

        # generate a list of files inside subpath
        flist = os.listdir(pp_subpath)

        # I want my flist to be just the tasks of interest. Other files like *.sh
        # files might exist in subpath, so here I need to find them and remove them
        # from flist before moving on. The remaining flist list is a list of strings
        # for each task directory
        pattern2find = "*.sh"
        matches = fnmatch.filter(flist, pattern2find)
        for match in matches:
            flist.remove(match)

        # loop over tasks in flist
        for task in flist:

            # start my bash script
            script = []
            script.append("#!/bin/bash")
            script.append("")
            script.append("# %s Plotting" % sub)
            script.append("")
            script.append("# %s" % task)

            # file to preprocess
            outpath = "%s/%s" % (pp_subpath, task)
            datafile1 = "%s/%s_%s%s.set" % (outpath,sub,task,fstem1)
            datafile2 = "%s/%s_%s%s.set" % (outpath,sub,task,fstem2)
            rawpath2use = "%s/%s" % (raw_subpath,task)
            rawfile = "%s/%s_%s.set" % (rawpath2use, sub, task)

            # MATLAB stuff to for running scorepochs
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_scorepochs('%s'); catch; cmi_scorepochs('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile1,datafile2)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

            # MATLAB stuff to make the plots
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_tsplot('%s'); catch; cmi_tsplot('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile1,datafile2)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_carpetplot('%s'); catch; cmi_carpetplot('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile1,datafile2)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_tsnrplot('%s'); catch; cmi_tsnrplot('%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,datafile1,datafile2)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

            # MATLAB stuff to compute data quality metrics
            matlab_command2run = "cd('%s'); maxNumCompThreads(%d); tic; try; cmi_preproc_quality('%s','%s'); catch; cmi_preproc_quality('%s','%s'); toc; exit; end; toc; exit;" % (codedir,maxNumCompThreads,rawfile,datafile1,rawfile,datafile2)
            script.append('matlab -nodesktop -r "%s"' % matlab_command2run)

            # Rmarkdown to make the report
            script.append("")
            r_code2run = "library(rmarkdown); setwd('%s'); render(input='%s/cmi_preproc_report.Rmd',output_file = '%s/%s_%s_preproc_report.html', params = list(subid = '%s', tasks = '%s'));" % (outpath, codedir,outpath,sub, task,sub, task)
            script.append('Rscript -e "%s"' % r_code2run)

            # Writing bash script to file ----------------------------------------------
            print("Saving preproc plots bash script for %s %s ..." %(sub,task))
            fname2save = "%s/%s/_02_preprocessing_plots_%s_%s.sh" % (pp_subpath,task,sub,task)
            fname = open(fname2save,"w")
            fname.write("\n".join(script)+"\n")
            fname.close()

            os.system("mkdir -p %s/_02_preprocessing_batches" % codedir)
            fname2save = "%s/_02_preprocessing_batches/_02_preprocessing_plots_%s_%s.sh" % (codedir,sub,task)
            fname = open(fname2save,"w")
            fname.write("\n".join(script)+"\n")
            fname.close()

            # if you're not just generating scripts and want to run the bash scripts
            # then do this below...
            if not script_only:

                if run_on_server:

                    print("Inserting batch into the queue ...")

                    # options for sbatch
                    # outlog = "%s/%s/%s_%s_preproc_out.log" % (pp_subpath,task,sub,task)
                    # errlog = "%s/%s/%s_%s_preproc_error.log" % (pp_subpath,task,sub,task)
                    sbatch_opts = "-J %s --tasks-per-node=1 --cpus-per-task=1 --time=02:00:00 --no-requeue --mem=8g" % (sub)

                    # run batch
                    os.system("sbatch %s %s" % (sbatch_opts,fname2save))

                else:

                    # run bash script
                    os.system("bash %s" % fname2save)
