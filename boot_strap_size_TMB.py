#coding=utf-8
import os,sys,random
import pandas as pd
import numpy as np
'''sys.argv[1] #只包含一列基因名称的文件；sys.argv[2] #每个基因梯度迭代次数;默认从100个gene开始迭代。梯度为50个gene'''
hg19_bed=os.getcwd()+'/cut_liftover_hg19_exon_current.gene.bed'  # hg19 bed
genepanel=sys.argv[1] #  一列基因
#dumplicate=sys.argv[1]
dup=1 
while dup<sys.argv[2]: #每个基因梯度迭代次数
	core_genelist_525=[]
	hg19_dic={}
	core_genelist2=[]
	'''
	with open(genepanel) as core_genelistfile:
		for eachline in core_genelistfile:
			if eachline.startswith('#'):
				continue
			each=eachline.strip().split('\t')[0]
			core_genelist2.append(each)
			if '/' not in each and each not in core_genelist_525:
				core_genelist_525.append(each)
			else:
				for eachx in each.split('/'):
					if eachx not in core_genelist_525:
						core_genelist_525.append(eachx)
	print(len(core_genelist_525),core_genelist_525[0])
	'''
	gene_mimic=[]
	with open(genepanel) as genepanelfile:
		for eachline in genepanelfile:
			gene_mimic.append(eachline.strip().split('\t')[0])
	WESgene=[]
	with open(hg19_bed) as hg19bed_file:
		for eachline in hg19bed_file:
			each=eachline.strip().split('\t')
			WESgene.append(each[-1])	
	for i in range(100,len(gene_mimic),50): #
		length=0;count=0;genex=[]
		core_genelist=[]
		addgene=random.sample(gene_mimic,i)
		for eachgene in addgene:
			if eachgene not in core_genelist:
				core_genelist.append(eachgene)
		print(len(core_genelist))
		outpanel=open('PANEL_dump_%s_%s'%(dup,len(core_genelist)),'w')
		for each in core_genelist:
			outpanel.write(each+'\n')
		outpanel.close()
		with open(hg19_bed) as hg19bed_file:
			for eachline in hg19bed_file:
				each=eachline.strip().split('\t')
				if each[-1] in core_genelist and each[-1] in  gene_mimic:
					length+=abs(int(each[2])-int(each[1]))+1
					if each[-1] not in genex:
						genex.append(each[-1])
		for each in core_genelist2:
			if '/' not in each:
				if each not in genex:
					print (each,)
			else:
				con=0
				for eachx in  each.split('/'):
					if eachx in genex:
						con+=1
				if con==0:
					print (each)
		if len(genex) !=len(core_genelist):
			print (len(genex),len(core_genelist))
			print ('注意panel中大约%s个gene不在hg19_bed文件中'%(len(core_genelist)-len(genex)))
			#exit('fuck dog!')		
		panel_exon=round(length/1000000.,3)
		print (panel_exon)
		###get SNV/INDEL in bed interval###
		WES_sample={}
		tumor_type=[]
		with open(os.getcwd()+'/ALL_8291_WES_panel_X_sample_TMB.xls') as WESfile:
			for eachline in WESfile:
				each=eachline.strip().split('\t')
				if each[1] not in tumor_type:
					tumor_type.append(each[1])
				name,types,WES_allcount,WES_count1,WES_TMB1,panel_ount1,panel_TMB1=[each[0],each[1],each[7],each[8],each[9],each[-4],each[-3]]
				WES_sample[name]='\t'.join([types,WES_allcount,WES_count1,WES_TMB1,panel_ount1,panel_TMB1])	
	
		varient_type1=['Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Missense_Mutation','Nonsense_Mutation']
		varient_type2=['Silent','Intron']
		out=open('WES_panel_TMB_dup_%s_%s'%(dup,len(core_genelist)),'w')
		out.write('sample\ttype\tnumber\tWES_count1\tWES_TMB1\t579_count\t579_TMB1\tpanelx_all\tpanelx_count\tpanelx_TMB\n')
		for eachsample in os.listdir(os.getcwd()+'/split_each_WES/'):
			sample=eachsample.split('_WES.maf')[0]
			mut_number=0;mut_numberx=0;mut_numbery=0
			mut_number_panel=0;mut_numberx_panel=0;mut_numbery_panel=0
			if sample in WES_sample.keys():
				with open(os.getcwd()+'/split_each_WES/'+eachsample) as MAFfile:
					for eachline in MAFfile:
						if 'Hugo_Symbol' in eachline:
							continue
						each=eachline.strip().split('\t')
						if each[0] in core_genelist:
							gene,chr,start,end,Variant_Classification,Tumor_Sample_Barcode,t_depth,t_ref_count,t_alt_count,CDS_position=[each[0],each[4],each[5],each[6],each[8],each[15],each[39],each[40],each[41],each[52]]
							freq=round(float(t_alt_count)/(int(t_ref_count)+int(t_alt_count)),3)
							mut_number+=1
							if Variant_Classification in varient_type1:
								mut_numberx+=1
							if Variant_Classification  in varient_type2:
								mut_numbery+=1
					WES_TMB1=round(mut_numberx/panel_exon,3);WES_TMB2=round(mut_numbery/panel_exon,3) #
					out.write(sample+'\t'+WES_sample[sample]+'\t'+str(mut_number)+'\t'+str(mut_numberx)+'\t'+str(WES_TMB1)+'\n')
		out.close()

		###=====cor==========
		os.system('mkdir %s'%('dup%s_panel_'%(dup)+str(len(core_genelist))))
		for eachtype in tumor_type:
			if 'type'  in eachtype:
				continue
			type_out=open('dup%s_panel_'%(dup)+str(len(core_genelist))+'/panel_'+str(len(core_genelist))+'_'+eachtype,'w')
			with open('WES_panel_TMB_dup_%s_%s'%(dup,len(core_genelist))) as panelx:
				for eachline in panelx:
					if eachline.startswith('sample'):
						type_out.write(eachline)
						continue
					if eachtype==eachline.strip().split('\t')[1]:
						type_out.write(eachline)
			type_out.close()
		cor_TMB=open('cor_TMB_dump_%s_%s'%(dup,str(len(core_genelist))),'w')
		cor_TMB.write('tumor type\tsperaman_TMB_x\tpearson_TMB_x\n')
		speraman_TMB={}
		for eachfile in os.listdir('dup%s_panel_'%(dup)+str(len(core_genelist))):
			data = pd.read_csv('dup%s_panel_'%(dup)+str(len(core_genelist))+'/'+eachfile,sep='\t')
			length=len(data)
			speraman_TMB_x=round(data['WES_TMB1'].corr(data['panelx_TMB'],'spearman'),3)
			speraman_TMB[eachfile]=speraman_TMB_x
			pearson_TMB_x=round(data['WES_TMB1'].corr(data['panelx_TMB'],'pearson'),3)
			cor_TMB.write(eachfile+'\t'+str(length)+'\t'+str(speraman_TMB_x)+'\t'+str(pearson_TMB_x)+'\n')	
		cor_TMB.close()
	dup+=1
	core_genelist=[]
