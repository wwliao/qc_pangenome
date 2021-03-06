# Quality Control of Pangenome Graph

## Preprocess pangenome VCF and merged VCF

1. Extract single-sample records from a multi-sample VCF

	```sh
	$ extract_single_sample_vcf.sh GRCh38-f1g-90-mc-aug11.vcf.gz HG00438
	```

2. Parse the single-sample VCF to get allele traversals

	```sh
	$ python3 get_allele_traversals.py HG00438.GRCh38-f1g-90-mc-aug11.vcf.gz
	```

3. Add stable allele traversals

	```sh
	$ python3 add_stable_traversals.py GRCh38-f1g-90-mc-aug11.gfa HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.txt
	```

4. Convert stable allele traversals to VCF format 

	```sh
	$ python3 convert_stable_traversals_to_vcf.py HG00438.header.txt GCA_000001405.15_GRCh38_no_alt_analysis_set.fa HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.txt
	```

5. Sort and remove duplicated variant records

	```sh
	$ bcftools sort -Ou HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.vcf | bcftools norm -d exact -Oz -o HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.sorted_rmdup.vcf.gz && bcftools index -t HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.sorted_rmdup.vcf.gz
	```

6. Convert symbolic VCF to explicit VCF

	```sh
	$ python3 convert_symbolic_to_explicit_vcf.py -s HG00438 GCA_000001405.15_GRCh38_no_alt_analysis_set.fa HG00438.merged.3.100.annot.final.reheader.vcf
	```

7. Sort and index merged variant records

	```sh
	$ bcftools sort -Oz -o HG00438.merged.3.100.annot.final.reheader.explicit.sorted.vcf.gz HG00438.merged.3.100.annot.final.reheader.explicit.vcf && bcftools index -t HG00438.merged.3.100.annot.final.reheader.explicit.sorted.vcf.gz
	```

## Evaluate performance for DEL >= 50bp

1. Extract DEL from VCF 

	```sh
	# Minigraph-Cactus VCF
	$ bcftools view -i "SVTYPE='DEL' & SVLEN<=-50" -Ou -o HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.sorted_rmdup.SV_DEL.vcf HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.sorted_rmdup.vcf.gz

	# Merged VCF
	$ bcftools view -i "SVTYPE='DEL' & SVLEN<=-50" -Ou -o HG00438.merged.3.100.annot.final.reheader.explicit.sorted.SV_DEL.vcf HG00438.merged.3.100.annot.final.reheader.explicit.sorted.vcf.gz
	```

2. Convert VCF to BED

	```sh
	$ python3 convert_vcf_to_bed.py HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.sorted_rmdup.SV_DEL.vcf
	$ python3 convert_merged_vcf_to_bed.py HG00438.merged.3.100.annot.final.reheader.explicit.sorted.SV_DEL.vcf
	$ cut -f 1-4 HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.stable.sorted_rmdup.SV_DEL.bed > mc.del.bed
	$ cut -f 1-4 HG00438.merged.3.100.annot.final.reheader.explicit.sorted.SV_DEL.bed > truth.del.bed
	```

3. Consider only variants in the dipcall confident regions

	```sh
	$ bedtools intersect -a mc.del.bed -b HG00438.f1_assembly_v2_genbank.dip.bed -u > mc.del.conf.bed
	$ bedtools intersect -a truth.del.bed -b HG00438.f1_assembly_v2_genbank.dip.bed -u > truth.del.conf.bed
	```

4. Consider only variants overlapping with at least 10% reciprocal overlap

	```sh
	$ bedtools intersect -a truth.del.conf.bed -b mc.del.conf.bed -f 0.1 -r -wa -wb > 10percent_roverlap.truth_query.bed
	$ cut -f 1-4 10percent_roverlap.truth_query.bed | sort -k1,1 -k2,2n -k3,3n -u > 10percent_roverlap.truth.bed
	$ cut -f 5-8 10percent_roverlap.truth_query.bed | sort -k1,1 -k2,2n -k3,3n -u > 10percent_roverlap.query.bed
	```

5. Compute coverage proportion 

	```sh
	$ bedtools coverage -a 10percent_roverlap.truth.bed -b 10percent_roverlap.query.bed > 10percent_roverlap.truth.cov.bed
	$ bedtools coverage -a 10percent_roverlap.query.bed -b 10percent_roverlap.truth.bed > 10percent_roverlap.query.cov.bed
	```

6. Calculate performance

	```sh
	$ python3 calc_performance.py -c 0.5 truth.del.bed mc.del.bed 10percent_roverlap.truth.cov.bed 10percent_roverlap.query.cov.bed
	```

# TODO

- Refine INV decomposition from stable allele traversals
- Solve COMPLEX decomposition from stable allele traversals
- Calculate SV calling performance metrics based on the GIAB stratification regions
- Manually check FP/FN variants using IGV

