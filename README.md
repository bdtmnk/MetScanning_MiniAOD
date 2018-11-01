# MetScanning
For recent instruction please visit: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETScanners
## Install
```
  cmsrel CMSSW_10_2_0
  cd CMSSW_10_2_0/src
  cmsenv
  git cms-init  
  git clone https://github.com/didukhle/MetScanning_MiniAOD.git
  scram b -j9
  ```
  You might need to run the following command if you want to access files via XROOT:
```
  voms-proxy-init --voms cms
```
## Run locally:
```
  cmsRun MetScanning/skim/python/skimMINIAOD.py
```

```
## Run with crab
In ``MetScanning/skim/crab/`` edit crab.py and adjust samples, JSON, and the EOS directory. 
Then do:
```
  cd MetScanning/skim/crab/
  python crab.py
```
In order to submit job with large input data:
```
  voms-proxy-init --voms cms
  source /cvmfs/cms.cern.ch/crab3/crab.sh
  
  Before submitting the jobs do a dryrun:

  crab --debug submit --config=data_crab.py --dryrun   

  In everything works fine you can then fully submit the jobs by:

  crab proceed
```
