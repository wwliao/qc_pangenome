# Quality Evaluation of Pangenome Graph

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
	$ add_stable_traversals.py GRCh38-f1g-90-mc-aug11.gfa HG00438.GRCh38-f1g-90-mc-aug11.allele_traversals.txt
	```