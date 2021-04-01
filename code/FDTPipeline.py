#! usr/bin/env python

# ## PIPELINE: FDTPipeline.py
# ## USAGE: python3 FDTPipeline --in=<inputs> --out=<outputs> [OPTIONS]
#    * requires python3 and FSL (calls FSL via python subprocess)
#
# ## Author(s)
#
# * Amy K. Hegarty, Intermountain Neuroimaging Consortium, University of Colorado Boulder
# * University of Colorado Boulder
#
# ## Product
#
# FSL Pipelines
#
# ## License
#
# <!-- References -->
# [FSL]: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki
# [pybids]: Yarkoni et al., (2019). PyBIDS: Python tools for BIDS datasets. Journal of Open Source Software, 4(40), 1294, https://doi.org/10.21105/joss.01294
#           Yarkoni, Tal, Markiewicz, Christopher J., de la Vega, Alejandro, Gorgolewski, Krzysztof J., Halchenko, Yaroslav O., Salo, Taylor, ? Blair, Ross. (2019, August 8). bids-standard/pybids: 0.9.3 (Version 0.9.3). Zenodo. http://doi.org/10.5281/zenodo.3363985
#

# ------------------------------------------------------------------------------
#  Show usage information for this script
# ------------------------------------------------------------------------------

def print_help():
  print("""
    Diffusion Preprocessing Pipeline
        Usage: """ + """ --in=<bids-inputs> --out=<outputs> [OPTIONS]
        OPTIONS
          --help                      show this usage information and exit
          --participant-label=        participant name for processing (pass only 1)
          --work-dir=                 (Default: /scratch) directory path for working 
                                        directory
          --clean-work-dir=           (Default: TRUE) clean working directory 
          --concat-preproc=           (Default: FALSE) boolean to select if input images
                                        should be concatenated before pre-processing
          --run-qc=                   (Default: TRUE) boolean to run automated quality 
                                        control for eddy corrected images
    ** fsl_sub used for job scheduleing. Set all necessary queue variables in system 
       environment before command execution 
        
        """)


# ------------------------------------------------------------------------------
#  Parse arguements for this script
# ------------------------------------------------------------------------------

def parse_arguments(argv):

    import os
    import sys
    import getopt

    #intialize arguements
    print("\nParsing User Inputs...")
    que = "ics"
    studyname = ""
    wd = "/working"
    cat = False
    qc = True
    cleandir = False


    try:
      opts, args = getopt.getopt(argv,"hi:o:",["in=","out=","help","participant-label=","work-dir=","clean-work-dir=","concat-preproc=","run-qc=","studyname="])
    except getopt.GetoptError:
      print_help()
      sys.exit(2)
    for opt, arg in opts:
      if opt in ("-h", "--help"):
         print_help()
         sys.exit()
      elif opt in ("-i", "--in"):
         inputs = arg
         if not os.path.exists(inputs):
           raise Exception("BIDS directory does not exist")
      elif opt in ("-o", "--out"):
         outputs = arg
      elif opt in ("--participant-label"):
         pid = arg
      elif opt in ("--work-dir"):
         wd = arg
      elif opt in ("--clean-work-dir"):
         cleandir = arg
      elif opt in ("--concat-preproc"):
         cat = arg
         if cat in ("TRUE", "True", "true"):
            cat = True
         elif cat in ("FALSE", "False", "false"):
            cat = False
         else:
            raise Exception("Error: --concat-preproc= [TRUE / FALSE]")
      elif opt in ("--run-qc"):
         qc = arg
         if qc in ("TRUE", "True", "true"):
            qc = True
         elif qc in ("FALSE", "False", "false"):
            qc = False
         else:
            raise Exception("Error: --run-qc= [TRUE / FALSE]")
      elif opt in ("--studyname"):
        studyname = arg
    if 'inputs' not in locals():
      print_help()
      raise Exception("Missing required argument --in=")
      sys.exit()
    if 'outputs' not in locals():
      print_help()
      raise Exception("Missing required argument --out=")
      sys.exit()
    if 'pid' not in locals():
      print_help()
      raise Exception("Missing required argument --participant-label=")
      sys.exit()
      
    print('Input Bids directory:\t', inputs)
    print('Derivatives path:\t', outputs)
    print('Participant:\t\t', str(pid))

    class args:
      def __init__(self, que, studyname, wd, inputs, outputs, pid, cat, qc, cleandir):
        self.que = que
        self.studyname = studyname
        self.wd = wd
        self.inputs = inputs
        self.outputs = outputs
        self.pid = pid
        self.concat= cat
        self.eddy_QC=qc
        self.cleandir=cleandir

    entry = args(que, studyname, wd, inputs, outputs, pid, cat, qc, cleandir)

    return entry

# ------------------------------------------------------------------------------
#  Parse Bids inputs for this script
# ------------------------------------------------------------------------------
def bids_data(entry):
    import os
    import glob
    import bids
    import json

    bids.config.set_option('extension_initial_dot', True)

    layout = bids.BIDSLayout(entry.inputs, derivatives=False, absolute_paths=True)

    if not os.path.exists(entry.outputs + '/FDT') or os.path.exists(entry.outputs + '/FDT/' + 'dataset_description.json'):
      os.makedirs(entry.outputs,exist_ok=True)
      os.makedirs(entry.outputs + '/FDT', exist_ok=True)

      # make dataset_description file...
      import json

      data = {
        'Name': 'FSL Diffusion Toolbox Minimal Preprocessing',
        "BIDSVersion": "1.1.1",
        "PipelineDescription": { 
              "Name": "FSL Diffusion Toolbox",
              "Version": "0.0.1",
              "CodeURL": "https://github.com/amyhegarty/FDT"
              },
        "CodeURL": "https://github.com/amyhegarty/FDT",
        "HowToAcknowledge": "Please cite all relevant works for FSL tools: topup, eddy, dtifit and python tools: pybids ( https://doi.org/10.21105/joss.01294,  https://doi.org/10.21105/joss.01294)"}

      with open(entry.outputs + '/FDT/' + 'dataset_description.json', 'w') as outfile:
          json.dump(data, outfile, indent=2)

    return layout

# ------------------------------------------------------------------------------
#  Main Pipeline Starts Here...
# ------------------------------------------------------------------------------

def worker(name,cmdfile):
    """Executes the bash script"""
    import subprocess
    from subprocess import PIPE
    process = subprocess.Popen(cmdfile.split(), stdout=PIPE, stderr=PIPE, universal_newlines=True)
    output, error = process.communicate()
    print(error)
    print('Worker: ' + name + ' finished')
    return

def run_topup(layout,entry):
    import os
    import sys
    import subprocess
    import multiprocessing
    
    # check if blip-up blip-down aquisition (ie dwi collected in opposing phase encoding directions)
    for dwi in layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'):
        ent = dwi.get_entities()
        if 'AP' in ent['direction']:
            if 'img1' not in locals():
                img1=dwi.path
                meta=dwi.get_metadata()
        elif 'PA' in ent['direction']:
            if 'img2' not in locals():
                img2=dwi.path

    # if blip-up blip-down is used: pull b0 from both images and merge for topup
    
    print('Applying topup: ')
    cmd='rm ' + entry.wd + '/acqparams.txt ; touch ' + entry.wd + '/acqparams.txt ;'
    if 'img1' in locals():
        cmd += 'fslroi ' + img1 + ' b0_AP 0 1; '  # takes the first volume of dwi image as b0
        cmd += 'echo "0 -1 0 ' + str(meta['TotalReadoutTime']) + '" >> ' + entry.wd + '/acqparams.txt ; ' # acqparameters for A -> P
        
        refimg = 'b0_AP'
        print('AP acquisition Using: ' + img1)

    if 'img2' in locals():
        cmd += 'fslroi ' + img2 + ' b0_PA 0 1; '  # takes the first volume of dwi image as b0  
        cmd += 'echo "0 1 0 ' + str(meta['TotalReadoutTime']) + '" >> ' + entry.wd + '/acqparams.txt ; ' # acqparameters for P -> A

        refimg = 'b0_AP'
        print('PA acquisition Using: ' + img2)

    if ('img1' in locals()) and ('img2' in locals()):
        cmd += 'fslmerge -t b0_APPA b0_AP b0_PA'
        refimg = 'b0_APPA'

    
    # write bash script for execution
    original_stdout = sys.stdout # Save a reference to the original standard output
    sys.stdout.flush()
    jid=[]

    with open(entry.wd + '/cmd_topup.sh', 'w') as fid:
      sys.stdout = fid # Change the standard output to the file we created.
      
      print('#!/usr/bin/bash')
      print('mkdir -p ' + entry.wd + '/topup')
      print('cd ' + entry.wd + '/topup')
      print(cmd)

      print("""topup --imain=""" + refimg + """ \
        --datain=../acqparams.txt \
        --config=b02b0.cnf \
        --out=topup_b0 \
        --iout=topup_b0_iout \
        --fout=topup_b0_fout
        --logout=topup""")

      sys.stdout = original_stdout # Reset the standard output to its original value

    # change permissions to make sure file is executable 
    os.chmod(entry.wd + '/cmd_topup.sh', 0o774)

    # run script
    cmd = "bash " + entry.wd + "/cmd_topup.sh"
    name = "topup"
    p = multiprocessing.Process(target=worker, args=(name,cmd))
    p.start()
    print(p)

    p.join()  # blocks further execution until job is finished

    ## end run_topup

# (Option 1) run eddy on each input scan seperately (no multi-scan concatination)
def run_eddy_opt1(layout,entry):
    import os
    import sys
    import subprocess
    import multiprocessing

    itr=0;
    nfiles = len(layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'))
    jobs=[];

    for dwi in layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'):
        if os.path.exists(entry.wd + '/eddy_dwi_' + str(itr) + '/eddy_unwarped_images.nii.gz'):
            itr=itr+1
            continue
            print(" ")

        img = dwi.path
        bval = layout.get_bval(dwi.path)
        bvec = layout.get_bvec(dwi.path)
        s=', '
        print('Using: ' + img)
        print('Using: ' + bval)
        print('Using: ' + bvec)

        # output filename...
        ent = layout.parse_file_entities(img)

        # Define the pattern to build out of the components passed in the dictionary
        pattern = "sub-{subject}/[ses-{session}/]sub-{subject}[_ses-{session}][_task-{task}][_acq-{acquisition}][_rec-{reconstruction}][_run-{run}][_echo-{echo}][_dir-{direction}][_space-{space}][_desc-{desc}]_{suffix}.nii.gz",

        # Add additional info to output file entities
        ent['space'] = 'native'
        ent['desc'] = 'preproc'

        outfile = layout.build_path(ent, pattern, validate=False, absolute_paths=False)
        outbval = outfile.replace('nii.gz','bval')
        outbvec = outfile.replace('nii.gz','bvec')
        outmask = outfile.replace('.nii.gz','_brain-mask.nii.gz')
        outqc   = outfile.replace('.nii.gz','.qc')

        print("corrected image: " + outfile)

        # write bash script for execution
        original_stdout = sys.stdout # Save a reference to the original standard output
        sys.stdout.flush()

        with open(entry.wd + '/cmd_eddy_dwi_' + str(itr) + '.sh', 'w') as fid:
          sys.stdout = fid # Change the standard output to the file we created.

          print('#!/usr/bin/bash')
          print('mkdir -p ' + entry.wd + '/eddy_dwi_' + str(itr))
          print('cd ' + entry.wd + '/eddy_dwi_' + str(itr))
          
          # add commands here...
          print('bvalfile=' + bval)
          print('n="$(grep -n "0" $bvalfile | cut -f1 -d":" )"')
          print('frame=`echo $n | cut -d" " -f1`')
          print('vol="$(($frame-1))"')
          print('echo "Using dwi volume : " $vol " for reference"')
          print('\n')
          print('fslroi ' + img + ' b0 $vol 1')
          topup_img = '../topup/topup_b0'
          acqparams = '../acqparams.txt'

          if 'AP' in img:
            inindex=1  # dwi images collected with acqparameters in row 1
          elif 'PA' in img: 
            inindex=2  # dwi images collected with acqparameters in row 1
          else:
            raise CustomError("Unable to determine if dwi image collected A->P or P->A")
          print('fslroi ' + img + ' b0 0 1') 
          print("""applytopup --imain=b0 \
                 --topup=""" + topup_img + """ \
                 --datain=""" + acqparams + """ \
                 --inindex=""" + str(inindex) + """ \
                 --method=jac \
                 --out=ref""")

          print('bet ref ref_brain -m -f 0.2')

          print('imglen=`fslval ' + img + ' dim4`')

          print('for c in $(seq 1 $imglen); do echo ' + str(inindex) + ' ; done > index.txt')

          print("""eddy_openmp --imain=""" + img + """ \
            --mask=ref_brain_mask \
            --index=index.txt \
            --acqp=""" + acqparams + """ \
            --bvecs=""" + bvec + """ \
            --bvals="""+ bval + """ \
            --fwhm=0 \
            --topup=""" + topup_img + """ \
            --flm=quadratic \
            --out=eddy_unwarped_images \
            --data_is_shelled""")

          if entry.eddy_QC == True:
            print("""eddy_quad eddy_unwarped_images  \
              -idx index.txt \
              -par """ + acqparams + """ \
              -m ref_brain_mask \
              -b """ + bval + """ \
              -g """ + bvec + """ \
              -f """ + topup_img + '_fout')

          print('mkdir -p $(dirname "' + entry.outputs + '/FDT/' + outfile + '")' )
          print('${FSLDIR}/bin/imcp eddy_unwarped_images.nii.gz ' + entry.outputs + '/FDT/' + outfile)
          print('cp -p ' + bval + ' ' + entry.outputs + '/FDT/' + outbval)
          print('cp -p ' + bvec + ' ' + entry.outputs + '/FDT/' + outbvec)
          print('cp -p ref_brain_mask.nii.gz ' + entry.outputs + '/FDT/' + outmask)
          print('cp -rp eddy_unwarped_images.qc/ ' + entry.outputs + '/FDT/' + outqc)

          print('ln -sf ' + bval + ' bvals')
          print('ln -sf ' + bvec + ' bvecs')
          print('ln -sf ' + img  + ' data.nii.gz')

          sys.stdout = original_stdout # Reset the standard output to its original value

        # change permissions to make sure file is executable 
        os.chmod(entry.wd + '/cmd_eddy_dwi_' + str(itr) + '.sh', 0o774)

        # run script
        cmd = 'bash ' + entry.wd + '/cmd_eddy_dwi_' + str(itr) + '.sh'
        name = 'eddy' + str(itr)
        p = multiprocessing.Process(target=worker, args=(name,cmd))
        jobs.append(p)
        p.start()

        itr = itr+1
        print(jobs)
    for job in jobs:
      job.join()  #wait for all eddy commands to finish

    ## end run_eddy_opt1

# (Option 2) run eddy on concatenated dwi scans (consistent to MRN, HCP pipelines (pairs))
def run_eddy_opt2(layout,entry):

    print('Concatenating dwi images...Running Eddy')

    import os
    import sys
    import subprocess
    import multiprocessing

    itr=0; s=', ';
    nfiles = len(layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'))

    # output filename...
    dwi = layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi')[0]
    ent = layout.parse_file_entities(dwi.path)

    # Define the pattern to build out of the components passed in the dictionary
    pattern = "sub-{subject}/[ses-{session}/]sub-{subject}[_ses-{session}][_task-{task}][_acq-{acquisition}][_rec-{reconstruction}][_run-{run}][_echo-{echo}][_dir-{direction}][_space-{space}][_desc-{desc}]_{suffix}.nii.gz",

    # Add additional info to output file entities
    ent['space'] = 'native'
    ent['desc'] = 'preproc'
    ent['direction'] = 'comb'  # combined set of dwi

    outfile = layout.build_path(ent, pattern, validate=False, absolute_paths=False)
    outbval = outfile.replace('nii.gz','bval')
    outbvec = outfile.replace('nii.gz','bvec')
    outmask = outfile.replace('.nii.gz','_brain-mask.nii.gz')
    outqc   = outfile.replace('.nii.gz','.qc')

    # get links to all input data...
    imglist=[]; bvallist=[]; bveclsit=[];

    # join bvals and bvecs from all scans
    cmd = ''
    cmd += 'mkdir -p ' + entry.wd + '/eddy_dwi_concat ; \n '
    cmd += 'cd ' + entry.wd + '/eddy_dwi_concat ; \n '
    cc=0
    for dwi in layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'):
      img = dwi.path
      print(img)
      cmd += 'ln -sf ' + dwi.path + ' dwi_' + str(cc) + '.nii.gz ; \n'
      cmd += 'ln -sf ' + layout.get_bval(dwi.path) + ' bval_' + str(cc) + ' ; \n' 
      cmd += 'ln -sf ' + layout.get_bvec(dwi.path) + ' bvec_' + str(cc) + ' ; \n' 

      topup_img = '../topup/topup_b0'
      acqparams = '../acqparams.txt'

      if 'AP' in img:
        inindex=1  # dwi images collected with acqparameters in row 1
      elif 'PA' in img: 
        inindex=2  # dwi images collected with acqparameters in row 1
      else:
        raise CustomError("Unable to determine if dwi image collected A->P or P->A")

      cmd += 'imglen=`fslval ' + img + ' dim4` ; \n'
      cmd += 'for c in $(seq 1 $imglen); do echo ' + str(inindex) + ' ; done > index_' + str(cc) + '.txt ; \n'
         
      cmd += 'fslroi ' + img + ' b0 0 1 ; \n'
      cmd += 'applytopup --imain=b0 --topup=' + topup_img + ' --datain=' + acqparams + ' --inindex=' + str(inindex) + ' --method=jac --out=ref ; \n'

      cmd += 'bet ref ref_brain -m -f 0.2 ; \n'
      cmd += 'mv ref_brain_mask.nii.gz ref_brain_mask_' + str(cc) + '.nii.gz ; \n'

      cc=cc+1

    cmd += 'touch bvecs; touch bvals; touch index_all.txt \n '
    cmd += 'paste -d " " $(ls bval_?) > bvals; \n '  # use to horizontally concatenate bval files
    cmd += 'paste -d " " $(ls bvec_?) > bvecs; \n '  # use to horizontally concatenate bvec files
    cmd += 'for i in $(ls index_?.txt); do cat $i ; done > index_all.txt; \n '

    # merge all raw images...
    cmd += 'images=$(ls dwi_?.nii.gz); \n '
    cmd += 'fslmerge -t data $images ; \n '

    # get mask...
    cmd += 'images=$(ls ref_brain_mask_?.nii.gz); \n '
    cmd += 'nscans=$(ls -1 ref_brain_mask_?.nii.gz | wc -l); \n'
    

    # do some fancy bash tricks to get a single line fslmaths command to get average mask
    cmd += 'echo "fslmaths ">tmp ; for i in $images; do echo "$i -add"; done >> tmp ; \n'  # print all individual masks to command
    cmd += "sed -i '$ s/.....$//' tmp; \n"                                                 # remove extra -add tag
    cmd += 'echo "-div $nscans avg_brain_mask" >> tmp; \n'                                 # add rest of command
    cmd += "a=`sed ':a;N;$!ba;s/" + "\\" + "n" + "/ /g' tmp`; \n"                          # make command single line
    cmd += "$a ; \n"                                                                       # run command
    # cmd += "cat tmp | tr " + repr('\n') + " ' ' > tmp; echo ' ' >> tmp; \n "
    # cmd += 'chmod a+x tmp; \n'
    # cmd += './tmp;  \n '
    cmd += 'fslmaths avg_brain_mask -thr 0.5 -bin brain_mask ; \n '

    # write bash script for execution
    original_stdout = sys.stdout # Save a reference to the original standard output
    sys.stdout.flush()
    jobs=[];

    with open(entry.wd + '/cmd_eddy_concat.sh', 'w') as fid:
      sys.stdout = fid # Change the standard output to the file we created.

      print('#!/usr/bin/bash')
      print(cmd)

      topup_img = '../topup/topup_b0'
      acqparams = '../acqparams.txt'

      print("""eddy_openmp --imain=data.nii.gz \
        --mask=brain_mask \
        --index=index_all.txt \
        --acqp=""" + acqparams + """ \
        --bvecs=bvecs \
        --bvals=bvals \
        --fwhm=0 \
        --topup=""" + topup_img + """ \
        --flm=quadratic \
        --out=eddy_unwarped_images \
        --data_is_shelled""")

      if entry.eddy_QC == True:
        print("""eddy_quad eddy_unwarped_images  \
          -idx index_all.txt \
          -par """ + acqparams + """ \
          -m brain_mask \
          -b bvals \
          -g bvecs \
          -f """ + topup_img + '_fout')

      print('mkdir -p $(dirname "' + entry.outputs + '/FDT/' + outfile + '")' )
      print('${FSLDIR}/bin/imcp eddy_unwarped_images.nii.gz ' + entry.outputs + '/FDT/' + outfile)
      print('cp -p bvals ' + entry.outputs + '/FDT/' + outbval)
      print('cp -p bvecs ' + entry.outputs + '/FDT/' + outbvec)
      print('cp -p brain_mask.nii.gz ' + entry.outputs + '/FDT/' + outmask)
      print('cp -rp eddy_unwarped_images.qc/ ' + entry.outputs + '/FDT/' + outqc)

      sys.stdout = original_stdout # Reset the standard output to its original value

    # change permissions to make sure file is executable 
    os.chmod(entry.wd + '/cmd_eddy_concat.sh', 0o774)

    # run script
    cmd = "bash " + entry.wd + "/cmd_eddy_concat.sh"
    name = 'eddy'
    p = multiprocessing.Process(target=worker, args=(name,cmd))
    jobs.append(p)
    p.start()

    print(jobs)

    for job in jobs:
      job.join()  #wait for all eddy commands to finish

    ## end run_eddy_opt2

# (Option 2) run tensor fitting on concatenated preprocessed image run with Opt 2
def run_dtifit_opt2(layout,entry):
    import os
    import sys
    import subprocess
    import multiprocessing

    # using the distortion corrected dti images, we will compute tensor values
    itr=0; s=', ';
    nfiles = len(layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'))
    jobs=[];

    for dwi in layout.get(subject=entry.pid, scope='derivatives', extension='nii.gz', direction='comb', suffix='dwi'):
        if os.path.exists(entry.wd + '/tensor_dwi_' + str(itr) + '/dwi_FA.nii.gz'):
          continue

        preproc_img = dwi.path
        bval = layout.get_bval(dwi.path)
        bvec = layout.get_bvec(dwi.path)
        mask = preproc_img.replace('.nii.gz','_brain-mask.nii.gz')

        print("Running dtifit: " + preproc_img)

        # write bash script for execution
        original_stdout = sys.stdout # Save a reference to the original standard output
        sys.stdout.flush()

        with open(entry.wd + '/cmd_tensor_dwi_' + str(itr) + '.sh', 'w') as fid:
          sys.stdout = fid # Change the standard output to the file we created.

          print('#!/usr/bin/bash')
          print('mkdir -p ' + entry.wd + '/tensor_dwi_' + str(itr))
          print('cd ' + entry.wd + '/tensor_dwi_' + str(itr))
          # add commands here...
          print("""dtifit --data=""" + preproc_img + """ \
            --mask=""" + mask + """ \
            --bvecs=""" + bvec + """ \
            --bvals="""+ bval + """ \
            --out=dwi """)

          #move outputs to derivative folder
          spath = preproc_img.replace("dwi.nii.gz","")
          print('for i in *.nii.gz; do ${FSLDIR}/bin/imcp $i ' + spath + '$i ; done')

          sys.stdout = original_stdout # Reset the standard output to its original value

        # change permissions to make sure file is executable 
        os.chmod(entry.wd + '/cmd_tensor_dwi_' + str(itr) + '.sh', 0o774)

        # run script
        cmd = 'bash ' + entry.wd + '/cmd_tensor_dwi_' + str(itr) + '.sh'
        name = 'dtifit' + str(itr)
        p = multiprocessing.Process(target=worker, args=(name,cmd))
        jobs.append(p)
        p.start()

        itr = itr+1
        print(jobs)
    for job in jobs:
      job.join()  #wait for all eddy commands to finish

    ## end run_dtifit_opt2

# (Option 1) concatenate pre-processed images (opt1) then run tensor fitting 
def run_dtifit_opt1(layout,entry):
    import os
    import sys
    import subprocess
    import multiprocessing

    # using the distortion corrected dti images, we will compute tensor values
    itr=0; s=', '
    nfiles = len(layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi'))
    jid=[]; jobs=[];

    # output filename...
    dwi = layout.get(subject=entry.pid, extension='nii.gz', suffix='dwi')[0]
    ent = layout.parse_file_entities(dwi.path)

    # Define the pattern to build out of the components passed in the dictionary
    pattern = "sub-{subject}/[ses-{session}/]sub-{subject}[_ses-{session}][_task-{task}][_acq-{acquisition}][_rec-{reconstruction}][_run-{run}][_echo-{echo}][_space-{space}][_desc-{desc}]_{suffix}.nii.gz",

    # Add additional info to output file entities
    ent['space'] = 'native'
    ent['desc'] = 'preproc'

    outfile = layout.build_path(ent, pattern, validate=False, absolute_paths=False)

    if not os.path.exists(entry.wd + '/tensor_dwi_join' + '/dwi_FA.nii.gz'):

      # join all eddy corrected images before calculating tensors...
      
      cmd = ''
      cmd += 'mkdir -p ' + entry.wd + '/tensor_dwi_join ; \n '
      cmd += 'cd ' + entry.wd + '/tensor_dwi_join ; \n '
      cmd += 'touch bvecs; touch bvals; touch index_all.txt \n '
      cmd += 'paste -d " " $(ls ../eddy_dwi_?/bvals) > bvals; \n '  # use to horizontally concatenate bval files
      cmd += 'paste -d " " $(ls ../eddy_dwi_?/bvecs) > bvecs; \n '  # use to horizontally concatenate bvec files
      cmd += 'for i in $(ls ../eddy_dwi_?/index.txt); do cat $i ; done > index_all.txt; \n '

      # merge all raw images...
      cmd += 'images=$(ls ../eddy_dwi_?/eddy_unwarped_images.nii.gz); \n '
      cmd += 'fslmerge -t data $images ; \n '

      # get average mask
      cmd += 'images=$(ls ../eddy_dwi_?/ref_brain_mask.nii.gz); \n '
      cmd += 'nscans=$(ls -1 ../eddy_dwi_?/ref_brain_mask.nii.gz | wc -l); \n'
      
      # do some fancy bash tricks to get a single line fslmaths command to get average mask
      cmd += 'echo "fslmaths ">tmp ; for i in $images; do echo "$i -add"; done >> tmp ; \n'  # print all individual masks to command
      cmd += "sed -i '$ s/.....$//' tmp; \n"                                                 # remove extra -add tag
      cmd += 'echo "-div $nscans avg_brain_mask" >> tmp; \n'                                 # add rest of command
      cmd += "a=`sed ':a;N;$!ba;s/" + "\\" + "n" + "/ /g' tmp`; \n"                           # make command single line
      cmd += "$a ; \n"                                                                       # run command
      # cmd += "cat tmp | tr " + repr('\n') + " ' ' > tmp; echo ' ' >> tmp; \n "
      # cmd += 'chmod a+x tmp; \n'
      # cmd += './tmp;  \n '
      cmd += 'fslmaths avg_brain_mask -thr 0.5 -bin brain_mask ; \n '
      
      print("Running dtifit: Joined DWI")

      # write bash script for execution
      original_stdout = sys.stdout # Save a reference to the original standard output
      sys.stdout.flush()

      with open(entry.wd + '/cmd_tensor_dwi_join.sh', 'w') as fid:
        sys.stdout = fid # Change the standard output to the file we created.

        print('#!/usr/bin/bash')
        
        # add commands here...
        print(cmd)

        print("""dtifit --data=data \
          --mask=brain_mask \
          --bvecs=bvecs \
          --bvals=bvals \
          --out=dwi""")

        #move outputs to derivative folder
        spath = outfile.replace("dwi.nii.gz","")
        print('for i in `ls *.nii.gz`; do ${FSLDIR}/bin/imcp $i ' + entry.outputs + '/FDT/' + spath + '$i ; done')

        sys.stdout = original_stdout # Reset the standard output to its original value

      # change permissions to make sure file is executable 
      os.chmod(entry.wd + '/cmd_tensor_dwi_join.sh', 0o774)

      # run script
      cmd = 'bash ' + entry.wd + '/cmd_tensor_dwi_join.sh'
      name = 'dtifit'
      p = multiprocessing.Process(target=worker, args=(name,cmd))
      jobs.append(p)
      p.start()

      print(jobs)

      for job in jobs:
        job.join()  #wait for all eddy commands to finish

    ## end run_dtifit_opt1

def run_cleanup(entry):

    import os
    import sys
    import subprocess
    import multiprocessing

    jobs=[];
    if entry.cleandir == True:
      # write bash script for execution
      original_stdout = sys.stdout # Save a reference to the original standard output
      sys.stdout.flush()

      with open(entry.wd + '/cmd_cleanup.sh', 'w') as fid:
        sys.stdout = fid # Change the standard output to the file we created.

        print('#!/usr/bin/bash')

        # add commands here...
        print('rm -Rf ' + entry.wd)

        sys.stdout = original_stdout # Reset the standard output to its original value

      # change permissions to make sure file is executable 
      os.chmod(entry.wd + '/cmd_cleanup.sh', 0o774)

      # run script
      cmd = 'bash ' + entry.wd + '/cmd_cleanup.sh'
      name = 'cleanup'
      p = multiprocessing.Process(target=worker, args=(name,cmd))
      jobs.append(p)
      p.start()

      print(jobs)

      for job in jobs:
        job.join()  #wait for all eddy commands to finish

    ## end run_cleanup


def main(argv):
    import glob
    import re
    import os
    import sys
    import warnings

    # get user entry
    entry = parse_arguments(argv)

    os.makedirs(entry.wd, exist_ok=True)
    logdir = entry.wd + '/logs'
    os.makedirs(logdir, exist_ok=True)

    # get participant bids path:
    bids = bids_data(entry)

    # pipeline: (1) topup, (2) eddy, (3) dtifit
    if not os.path.exists(entry.wd + '/topup/topup_b0_iout.nii.gz'):
        run_topup(bids,entry)
    
    # two run options: 
    if entry.concat == False:
      # (1) distortion and eddy correct each aquisition seperately
      run_eddy_opt1(bids,entry)

      bids.add_derivatives(entry.outputs + '/FDT')

      run_dtifit_opt1(bids,entry)

    else:
      # (2) concatenate all aquisitions before preprocessing
      run_eddy_opt2(bids,entry)

      bids.add_derivatives(entry.outputs + '/FDT')

      run_dtifit_opt2(bids,entry)

    # clean-up
    run_cleanup(entry)
    

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
