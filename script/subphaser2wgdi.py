#coding: utf-8
import sys
import pandas as pd

#print("sys.argv[1]形式例子：1A	0   60000   C")

data = pd.read_csv(sys.argv[1], sep="\t", header=None)
data.columns = ['chr','start','end','exchangefr']

sgs = set(data['exchangefr'])
d_sgs = {sg: i+1 for i, sg in enumerate(sorted(sgs))}

datagff = pd.read_csv(sys.argv[2], sep="\t", header=None)
datagff.columns = ['chr','gene','start','end','chainorder','geneorder','gene']
datagff['chr'] = datagff['chr'].astype(dtype='str')
databychr = datagff.groupby('chr')

#print(type(databychr.get_group('14').reset_index()))

#print(databychr.size())

'''
chr14=pd.DataFrame(databychr.get_group('14'))
for i in chr14.index:
	print(i)
'''
#print(data.loc[0,'exchangefr'])
out=pd.DataFrame()

for i in range(0,len(data)):
	start = data['start'][i] 
	end = data['end'][i]
	chrid = str(data['chr'][i]) #染色体id
	datachrid=databychr.get_group(chrid) #根据染色体id分组

#	print(datachrid['start'].tolist()[0:3])
	locget=datachrid['start'].tolist() #将该染色体gff的基因start列加入bed文件范围左侧做成列表
#	print(type(locget))
	locget.append(start)
#	print(sorted(locget).index(start)) #得到bed文件范围左侧在列表位置
	locbam=sorted(locget).index(start)
	if (locbam==0):
		geneorderleft=datachrid['geneorder'][datachrid.index[0]+locbam]
	elif (datachrid['end'][datachrid.index[0]+locbam-1]>start): #判断bed文件碱基范围左侧的上个基因的右侧end是否超过bed文件碱基范围左侧位置
		#print(datachrid['geneorder'][datachrid.index[0]+locbam-1]) 
		geneorderleft=datachrid['geneorder'][datachrid.index[0]+locbam-1]
	elif (len(datachrid)<=locbam): 
		#print(datachrid['geneorder'][datachrid.index[0]+locbam-1])
		geneorderleft=datachrid['geneorder'][datachrid.index[0]+locbam-1]
	else:
		#print(datachrid['geneorder'][datachrid.index[0]+locbam])
		geneorderleft=datachrid['geneorder'][datachrid.index[0]+locbam]

	locget2=datachrid['end'].tolist() 
	locget2.append(end)
	locbam2=sorted(locget2).index(end)
	if (len(datachrid)<=locbam2): #判断bed右侧结束位置是否是该染色体注释最后一行，如果这样，下一行，即datachrid.index[0]+locbam2不存在
		#print(datachrid['geneorder'][datachrid.index[0]+locbam2-1])
		geneorderright=datachrid['geneorder'][datachrid.index[0]+locbam2-1]
	elif (datachrid['start'][datachrid.index[0]+locbam2]<end): #判断bed文件碱基范围右侧的下个基因的左侧start是否小于bed文件碱基范围右侧位置
		#print(datachrid['geneorder'][datachrid.index[0]+locbam2])
		geneorderright=datachrid['geneorder'][datachrid.index[0]+locbam2]
	else:
		#print(datachrid['geneorder'][datachrid.index[0]+locbam2-1])
		geneorderright=datachrid['geneorder'][datachrid.index[0]+locbam2-1]
#	print(type(geneorderleft))
	#print(str(data.loc[i,'chr'])+"\t"+str(geneorderleft)+"\t"+str(geneorderright)+"\t"+str(data.loc[i,'exchangefr']))
	sg = str(data.loc[i,'exchangefr'])
	sgi = d_sgs[sg]
	outdf=pd.DataFrame([ [str(data.loc[i,'chr']),str(geneorderleft),str(geneorderright),'red', sgi] ],columns=['chr','start','end','color','sg'])
#	print(outdf.loc[0,:])
	out=out.append(outdf,ignore_index=True)

#print(out.head(5))
out.to_csv(sys.argv[3],sep="\t", index=False,header=None) 
