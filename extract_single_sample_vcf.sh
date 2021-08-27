#!/bin/bash
VCF=$1
SAMPLE=$2
PREFIX=$(basename $VCF .vcf.gz)

bcftools view -a -s $SAMPLE -Ou $VCF \
    | bcftools view -e 'GT="ref" | GT=".|."' -Oz -o ${SAMPLE}.${PREFIX}.vcf.gz \
	&& bcftools index -t ${SAMPLE}.${PREFIX}.vcf.gz
