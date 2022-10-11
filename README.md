# HeavyFlavorAnalysis

This repository is about the merging of the SKIM and Ntuple steps for the double Jpsi cross section measurement study.  
The SKIM is from:  
https://github.com/zhenhu/HeavyFlavorAnalysis  
The core file of the Ntuple Maker is: HeavyFlavorAnalysis/Onia2MuMu/src/MuMuGammaRootupler.cc:  
https://github.com/JinfengLiu97/HeavyFlavorAnalysis/blob/main/Onia2MuMu/src/MuMuGammaRootupler.cc  
****
    cmsrel CMSSW_10_6_20
    cd CMSSW_10_6_20/src/
    cmsenv
    git clone https://github.com/JinfengLiu97/HeavyFlavorAnalysis
    scram b
    cd HeavyFlavorAnalysis/Onia2MuMu/test/
    voms-proxy-init -voms cms --valid 172:00  
Do the local test first:  
`cmsRun BPH_SKIM_Ntuple_data.py`  
After the test is passed, submit the CRAB request by:  
`crab submit crab3_SKIM_Ntuple.py`  
All the datasets and JSONs are already listed in the CRAB configuration file, modify them properly before the submission.  

## About the Ntuple

As for now, the Ntuple Maker mainly saves four variables about the distance (prompt/non-prompt distinguishment).  
They and corresponding branches are listed here:  
* ctau: fourMuFit_ups1/2_cTau_MC/noMC  
* LxyPV: fourMuFit_ups1/2_LxyPV_MC/noMC  
* Significance of LxyPV: fourMuFit_ups1/2_LxyPVSig_MC/noMC  
* Significance of the distance of Jpsis decay vertexes (dJ/psi): fourMuFit_DistanceSig_MC/noMC  

You can also find related variables according to these branches.  

## Wish you a nice journey in data processing!

