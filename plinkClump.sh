

module load plink

# clump AEI
plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/geuvadisAEIgwas.mlma_e-05_Hits.txt --clump-p1 5e-05 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/geuvadisAEIgwas.mlma_e-05_clumpedHits.txt

plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/geuvadisAEIgwas.mlma_e-08_Hits.txt --clump-p1 5e-08 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/geuvadisAEIgwas.mlma_e-08_clumpedHits.txt

# clump CEI
plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/geuvadisCEIgwas.mlma_e-05_Hits.txt --clump-p1 5e-05 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/geuvadisCEIgwas.mlma_e-05_clumpedHits.txt

plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/geuvadisCEIgwas.mlma_e-08_Hits.txt --clump-p1 5e-08 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/geuvadisCEIgwas.mlma_e-08_clumpedHits.txt

# clump MEL
plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/meanEditLevelGWAS.mlma_e-05_Hits.txt --clump-p1 5e-05 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/meanEditLevelGWAS.mlma_e-05_clumpedHits.txt

plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/meanEditLevelGWAS.mlma_e-08_Hits.txt --clump-p1 5e-08 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/meanEditLevelGWAS.mlma_e-08_clumpedHits.txt

# clump McEL
plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/cdsEditLevelGWAS.mlma_e-05_Hits.txt --clump-p1 5e-05 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/cdsEditLevelGWAS.mlma_e-05_clumpedHits.txt

plink --bfile genoFiles/allGeuvadisSampGeno \
--clump assocRes/cdsEditLevelGWAS.mlma_e-08_Hits.txt --clump-p1 5e-08 --clump-p2 0.01 \
--clump-r2 0.1 --clump-kb 250 \
--out assocRes/cdsEditLevelGWAS.mlma_e-08_clumpedHits.txt
