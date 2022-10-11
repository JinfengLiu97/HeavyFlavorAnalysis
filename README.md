# HeavyFlavorAnalysis

This repository is about the merging of the SKIM and Ntuple steps for the double Jpsi cross section measurement study.  
The SKIM is from: https://github.com/zhenhu/HeavyFlavorAnalysis  
The core file of the Ntuple Maker is: HeavyFlavorAnalysis/Onia2MuMu/src/MuMuGammaRootupler.cc  
****
    cmsrel CMSSW_10_6_20
    cd CMSSW_10_6_20/src/
    cmsenv
    git clone ???
    scram b
    cd HeavyFlavorAnalysis
    cd Onia2MuMu/test/
    voms-proxy-init -voms cms --valid 172:00  
Do the local test first:  
`cmsRun BPH_SKIM_Ntuple_data.py`  
After the test is passed, submit the CRAB request by:  
`crab submit crab3_SKIM_Ntuple.py`  
All the datasets and JSONs are already listed in the CRAB configuration file, modify them properly before the submission.  
Wish you a nice journey in data processing!
