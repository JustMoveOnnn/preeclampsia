java -cp gsea-3.0.jar -Xmx4096m xtools.gsea.Gsea -out deci_Endo/ -rnd_seed 13579 -set_min 5 -set_max 1000 -res deci_Endo_EXP.txt -cls deci_Endo_PT.cls -gmx c2.all.v7.0.symbols.gmt -plot_top_x 2000 -permute gene_set -collapse false