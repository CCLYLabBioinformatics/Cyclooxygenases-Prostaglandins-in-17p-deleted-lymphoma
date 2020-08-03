######------GSEA_shAlox15bVshREN_vehicl#####
java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/h.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_h -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c1.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c1 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c2.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c2 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c3.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c3 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c4.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c4 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c5.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c5 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c6.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c6 -gui false

java -cp /mnt/data/userdata/xiangyu/programme/gsea/gsea-3.0.jar -Xmx1024m xtools.gsea.Gsea \
-res /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/Alox.gct \
-cls /mnt/data/userdata/cclabshare/QL/sortedByCoord.out.bam/GSEA/samplename.cls#shAlox15b_versus_shREN \
-gmx /mnt/data/userdata/xiangyu/programme/gsea/msigdb_v6.1_files_to_download_locally/msigdb_v6.1_GMTs/c7.all.v6.1.symbols.gmt \
-chip /mnt/data/userdata/xiangyu/programme/gsea/chip/ENSEMBL_mouse_gene.chip \
-collapse true -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set \
-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis \
-metric Signal2Noise -sort real -order descending -include_only_symbols true \
-make_sets true -median false -num 100 -plot_top_x 50 -rnd_seed timestamp \
-save_rnd_lists false -set_max 500 -set_min 15 -zip_report false \
-out ./shAlox15b_versus_shREN_c7 -gui false


