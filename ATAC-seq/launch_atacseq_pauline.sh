nextflow run nf-core/atacseq \
-r 1.2.1 \
-resume atac_seq_3 \
-profile docker \
-work-dir /media/legrand-lab/Stockage/InstitutNeuroMyoGene_DATA/EqLegrand/atacseq_pauline/work \
-params-file nf-params_atacseq.json
