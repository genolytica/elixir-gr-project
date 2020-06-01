## Basic example

```
Rscript run_metaseqR2.R \
  --targets=my_targets.txt \
  --contrast=A_vs_B \
  --org=hg19 \
  --counttype=exon \
  --normalization=edger \
  --statistics=edger \
  --figformat=png \
  --xprtwhere=. \
  --rc=0.5 \
  --exportr2c
```
