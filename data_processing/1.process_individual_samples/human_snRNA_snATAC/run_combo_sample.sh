eval "$(conda shell.bash hook)"
export TMPDIR=/diskmnt/Projects/SenNet_primary/tmp/
sample=$1
cd /diskmnt/Projects/SenNet_primary/Cellranger-arc/
if [ -d "$sample"/outs ]
then
echo echo "Cellranger for $sample exists, moving to next step"
else
cellranger-arc count \
--id ${sample} \
--reference /diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--libraries /diskmnt/Projects/SenNet_primary/Cellranger-arc/preprocessing/${sample}/library.csv \
--disable-ui --localcores 40 --localmem 100
fi

cd /diskmnt/Projects/SenNet_primary/Scrublet
if [ -f combo_merged/${sample}/${sample}_combo_scrublet_output_table.csv ]
        then
                echo "Scrublet for $sample exists, moving to next step"
        else

cd /diskmnt/Projects/SenNet_primary/Scrublet/RNA
conda activate scrublet&&bash /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/Scrublet/scripts/scrublet-RNA-one-sample.sh \
/diskmnt/Projects/SenNet_primary/Cellranger-arc/ $sample

cd /diskmnt/Projects/SenNet_primary/Scrublet/ATAC
bash /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/Scrublet/scripts/scrublet-ATAC-one-sample.sh \
/diskmnt/Projects/SenNet_primary/Cellranger-arc/ $sample

cd /diskmnt/Projects/SenNet_primary/Scrublet/combo_merged
bash /diskmnt/Projects/Users/austins2/tools/scrublet-auto-combining.sh \
/diskmnt/Projects/SenNet_primary/Cellranger-arc/ \
/diskmnt/Projects/SenNet_primary/Scrublet/RNA/ \
/diskmnt/Projects/SenNet_primary/Scrublet/ATAC/
fi

cd /diskmnt/Projects/SenNet_primary/rds_objects/
conda activate signac1.8&&Rscript /diskmnt/Projects/Users/allakarpova/scripts/multiome_sc_analysis/seurat_pipeline_v5.2.combo.R \
-s ${sample} \
-d /diskmnt/Projects/SenNet_primary/Cellranger-arc/${sample} \
-m /diskmnt/Projects/Users/allakarpova/Tools/miniconda3/envs/signac1.8/bin/macs2 \
-o /diskmnt/Projects/SenNet_primary/rds_objects/ \
--chrom_size /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt \
--prf_min 1000 --pct_min 15 --ns_max 5 --pc_first 2 --pc_num 50 \
--scrublet /diskmnt/Projects/SenNet_primary/Scrublet/combo_merged/${sample}/

