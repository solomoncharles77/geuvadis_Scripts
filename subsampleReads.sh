

zcat A1750_1.fq.gz | awk 'NR<=4000000' | gzip > A1750_1_sub.fq.gz
zcat A1750_2.fq.gz | awk 'NR<=4000000' | gzip > A1750_2_sub.fq.gz

zcat A1802_1.fq.gz | awk 'NR<=4000000' | gzip > A1802_1_sub.fq.gz
zcat A1802_2.fq.gz | awk 'NR<=4000000' | gzip > A1802_2_sub.fq.gz

zcat A1807_1.fq.gz | awk 'NR<=4000000' | gzip > A1807_1_sub.fq.gz
zcat A1807_2.fq.gz | awk 'NR<=4000000' | gzip > A1807_2_sub.fq.gz
