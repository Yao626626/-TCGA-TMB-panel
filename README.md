# -TCGA-TMB-panel
利用TCGA数据迭代设计TMB检测panel
使用： python boot_strap_size_TMB.py  {MSK_FMI_525.txt}  {1021}
默认已写入脚本中的文件：
1. hg19的bed文件：cut_liftover_hg19_exon_current.gene.bed；
2. 32个癌种，8291个TCGA样本的TMB总文件：ALL_8291_WES_panel_X_sample_TMB.xls
3. WES数据高可信度结果maf文件 百度网盘链接：https://pan.baidu.com/s/1O0EUDZodZ0sT-8Qc8MQDDA  提取码：wsxt 

其它输入文件：仅含一列基因名称的文件，如 MSK_FMI_525.txt，并指定迭代次数：如 1021

输出结果：TMB相关分析结果文件：cor_TMB_dump_{迭代次数}_{基因数目}；PANEL基因列表文件：PANEL_dump_99_625；原始TMB信息汇总文件：WES_panel_TMB_dup_99_982；32个癌种中样本TMB信息文件夹：dup9_panel_974
