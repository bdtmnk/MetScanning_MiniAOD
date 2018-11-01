# MetScanning
For recent instruction please visit: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETScanners
## Install
```
  cmsrel CMSSW_10_1_2_patch2
  cd CMSSW_10_1_2_patch2/src
  cmsenv
  git cms-init
  git cms-addpkg RecoMET/METFilters 
  
  git clone https://github.com/didukhle/MetScanning.git
  scram b -j9
  ```
  You might need to run the following command if you want to access files via XROOT:
```
  voms-proxy-init --voms cms
```
## Run on local file:
```
  cmsRun MetScanning/skim/python/skim.py
```
## Run on AOD files:
```
  cmsRun MetScanning/skim/python/skimAOD.py
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

  crab --debug submit --config=crab.py --dryrun   

  In everything works fine you can then fully submit the jobs by:

  crab proceed
```
