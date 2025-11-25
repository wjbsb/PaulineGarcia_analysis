#mkdir -p "/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/FIXED"

#bash /home/mna_bioinfo/Bureau/installation/BBMap_39.06/bbmap/repair.sh \
#-in1="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/Input_WT1_S13_L002_R1_001.fastq.gz" \
#-in2="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/Input_WT1_S13_L002_R2_001.fastq.gz" \
#-out1="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/FIXED/Input_WT1_S13_L002_R1_001.fastq.gz" \
#-out2="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/FIXED/Input_WT1_S13_L002_R2_001.fastq.gz" \
#-outs="/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/FIXED/singletons_Input_WT1_S13_L002_R2_001.fastq.gz" \ 
#repair

for fastqgz in `find /media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated -name "*R1_001.fastq.gz"`; do
	
	f=`basename $fastqgz`
	
	r=`cut -d "_" -f5 <<< $f`
	name=`cut -d "_" -f1-4 <<< $f`
	
	r1=/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/${name}_${r}_001.fastq.gz
	r2=/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FastQLanesSeparated/${name}_R2_001.fastq.gz
	
	o1=/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FIXED/${name}_${r}_001.fastq.gz
	o2=/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FIXED/${name}_R2_001.fastq.gz
	s=/media/mna_bioinfo/MNA2_Stockage/FABIEN_LEGRAND/PAULINE_GARCIA/chipseq/transfer_6412495_files_71e4e46a/FIXED/singletons_${name}.fastq.gz
	
	bash /home/mna_bioinfo/Bureau/installation/BBMap_39.06/bbmap/repair.sh \
	-in1=${r1} \
	-in2=${r2} \
	-out1=${o1} \
	-out2=${o2} \
	-outs=${s} \ 
	repair
done
