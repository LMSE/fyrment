import cobra
import pandas as pd
from libsbml import *

root_aybrah='../aybrah/'

# local_path
#root_aybrah='https://github.com/kcorreia/aybrah/raw/master/'


path_aybrah_xlsx=root_aybrah+'aybrah.xlsx'


taxonomy=pd.read_excel(path_aybrah_xlsx,sheet_name='taxon_nodes').fillna('')

taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))[1]

taxonomy['genes']=''
taxonomy['reactions']=''
taxonomy['metabolites']=''

for index,row in taxonomy.iterrows():
	name=taxonomy['name'][index]
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	readable=taxonomy['readable'][index].replace('.','_')
	path_xml_dir = './genre/xml/'+str(level_order).zfill(2)+'_'+level.lower()
	path_xml_file=path_xml_dir+'/i'+readable+'.xml'
	oid=taxonomy['oids'][index]
	if oid in ['cpr','ani']:
		continue
	model=cobra.io.read_sbml_model(path_xml_file)
	reactions=model.reactions
	metabolites=model.metabolites
	genes=model.genes
	#
	taxonomy['genes'][index]=len(genes)
	taxonomy['reactions'][index]=len(reactions)
	taxonomy['metabolites'][index]=len(metabolites)


taxonomy.to_csv('fyrment_stats.txt',sep='\t',index=False)



for index,row in taxonomy[taxonomy.level=='Strain'].iterrows():
	#break
	oid=taxonomy['oids'][index]
	name=' '.join(taxonomy['name'][index].split(' ')[:2])
	row=oid+' | '+name+'\t'+','.join([str(taxonomy[col][index]) for col in 'genes|reactions|metabolites'.split('|')])
	print(row)



	