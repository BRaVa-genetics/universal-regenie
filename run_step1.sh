bash 01_step1_fitNULLGLMM.sh -t binary --plink ukb_array.wes_450k_qc_pass_eur.for_vr.chr21 --sparseGRM ukb_array.wes_450k_qc_pass_allpop.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx --sparseGRMID ukb_array.wes_450k_qc_pass_allpop.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt --phenoFile UKBB_WES_Infertility_EUR_ONLY.txt --phenoCol "female_infertility_binary" --covarColList "age,assessment_centre" --categCovarColList "assessment_centre" --sampleIDCol "IID" --outputPrefix "out/" --isSingularity true 
