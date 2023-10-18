#!/bin/tcsh

fastq-dump --split-3 SRR4733580 > /dev/null
echo "SRR4733580 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733580_1.fastq  SRR4733580.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733580.bam SRR4733580.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733581 > /dev/null
echo "SRR4733581 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733581_1.fastq  SRR4733581.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733581.bam SRR4733581.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733582 > /dev/null
echo "SRR4733582 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733582_1.fastq  SRR4733582.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733582.bam SRR4733582.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733583 > /dev/null
echo "SRR4733583 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733583_1.fastq  SRR4733583.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733583.bam SRR4733583.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733584 > /dev/null
echo "SRR4733584 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733584_1.fastq  SRR4733584.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733584.bam SRR4733584.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733585 > /dev/null
echo "SRR4733585 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733585_1.fastq  SRR4733585.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733585.bam SRR4733585.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733586 > /dev/null
echo "SRR4733586 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733586_1.fastq  SRR4733586.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733586.bam SRR4733586.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733587 > /dev/null
echo "SRR4733587 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733587_1.fastq  SRR4733587.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733587.bam SRR4733587.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733588 > /dev/null
echo "SRR4733588 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733588_1.fastq  SRR4733588.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733588.bam SRR4733588.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733589 > /dev/null
echo "SRR4733589 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733589_1.fastq  SRR4733589.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733589.bam SRR4733589.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733590 > /dev/null
echo "SRR4733590 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733590_1.fastq  SRR4733590.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733590.bam SRR4733590.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733591 > /dev/null
echo "SRR4733591 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733591_1.fastq  SRR4733591.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733591.bam SRR4733591.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733592 > /dev/null
echo "SRR4733592 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733592_1.fastq  SRR4733592.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733592.bam SRR4733592.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733593 > /dev/null
echo "SRR4733593 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733593_1.fastq  SRR4733593.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733593.bam SRR4733593.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733594 > /dev/null
echo "SRR4733594 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733594_1.fastq  SRR4733594.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733594.bam SRR4733594.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733595 > /dev/null
echo "SRR4733595 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733595_1.fastq  SRR4733595.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733595.bam SRR4733595.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733596 > /dev/null
echo "SRR4733596 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733596_1.fastq  SRR4733596.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733596.bam SRR4733596.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733597 > /dev/null
echo "SRR4733597 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733597_1.fastq  SRR4733597.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733597.bam SRR4733597.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733598 > /dev/null
echo "SRR4733598 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733598_1.fastq  SRR4733598.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733598.bam SRR4733598.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733599 > /dev/null
echo "SRR4733599 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733599_1.fastq  SRR4733599.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733599.bam SRR4733599.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733600 > /dev/null
echo "SRR4733600 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733600_1.fastq  SRR4733600.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733600.bam SRR4733600.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733601 > /dev/null
echo "SRR4733601 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733601_1.fastq  SRR4733601.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733601.bam SRR4733601.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

fastq-dump --split-3 SRR4733602 > /dev/null
echo "SRR4733602 (paired)" >> bowtie_test.output
bowtie -p10 --best --strata --no-unal -m 1 -X 50000 -5 10  -3 20 -n 3 --sam --chunkmbs 1000 /ReplaceRootPath/Databases/Drosophila/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex/genome SRR4733602_1.fastq  SRR4733602.sam |& tee -a bowtie.output

samtools view -bS -o BAM/SRR4733602.bam SRR4733602.sam > /dev/null

rm SRR*.sam
rm SRR*.fastq

