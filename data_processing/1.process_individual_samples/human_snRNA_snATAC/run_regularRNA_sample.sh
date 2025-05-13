eval "$(conda shell.bash hook)"
export TMPDIR=/diskmnt/Projects/SenNet_primary/tmp/
sample=$1
cd /diskmnt/Projects/SenNet_primary/Cellranger/
if [ -d "$sample"/outs ]
        then
		echo "Cellranger for $sample exists, moving to next step"
else
	cellranger count \
	--id ${sample} \
	--transcriptome /diskmnt/Datasets/Reference/Cellranger-2020-A/refdata-gex-GRCh38-2020-A \
	--fastqs /diskmnt/Projects/SenNet_primary/Cellranger/preprocessing/${sample}/ \
	--localcores 20 --localmem 100 --include-introns
fi

cd /diskmnt/Projects/SenNet_primary/Scrublet/RNA
if [ -d "$sample" ]
	then
		echo "Scrublet for $sample exists, moving to next step"
	else
		conda activate scrublet&&bash /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/Scrublet/scripts/scrublet-RNA-one-sample.sh /diskmnt/Projects/SenNet_primary/Cellranger/ $sample
fi

cd /diskmnt/Projects/SenNet_primary/rds_objects/
conda activate signac1.8&&Rscript /diskmnt/Projects/Users/allakarpova/scripts/snRNA_my/Seurat_v4_analysis_auto_v1.3.R \
-i /diskmnt/Projects/SenNet_primary/Cellranger/${sample}/outs/raw_feature_bc_matrix/ \
-s ${sample} -o /diskmnt/Projects/SenNet_primary/rds_objects/ \
--scrublet /diskmnt/Projects/SenNet_primary/Scrublet/RNA/${sample}/

