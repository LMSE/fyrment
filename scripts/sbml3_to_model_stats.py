
import pandas as pd
import xml.etree.ElementTree as ET

taxonomy=pd.read_csv('taxonomy_level_oids.txt',sep='\t')


drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)

taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))[1]


organisms=pd.read_excel('/Volumes/5TB/Organized/github/aybrah_old/aybrah_0.1.8_production.xlsx',sheet_name='curated_taxonomy_fungi').fillna('').set_index('order')



version_aybraham='0.0.7'
model_xlsx='aybraham_v0.0.7.xlsx'

rxns=pd.read_excel(model_xlsx,sheet_name='reactions',skiprows=[0]).fillna('').set_index('Rxn name')

subsystems=sorted(set([(subsystem.split('|')[0],subsystem.split('|')[1]) for subsystem in rxns['Subsystem'].tolist()]))
output=[]
for order,s in subsystems:
	if s not in output:
		output.append(s)

subsystems=output




model_stats=pd.DataFrame(columns='genes|reactions|metabolites'.split('|'))

model_stats_subsystem_rxns=pd.DataFrame(columns=subsystems)
model_stats_subsystem_genes=pd.DataFrame(columns=subsystems)

len(model_stats_subsystem_rxns.columns)

for index,row in taxonomy.iterrows():
#for index,row in taxonomy[taxonomy['level']=='Strain'].iterrows():
	#
	name=taxonomy['name'][index]
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	readable=taxonomy['readable'][index]
	if level=='Strain':
		print name
	else:
		continue
	oid=taxonomy['oids'][index]
	path_xml_dir = './xml/'+str(level_order).zfill(2)+'_'+level.lower()
	path_xml_file=path_xml_dir+'/i'+readable+'.xml'
	#
	tree = ET.parse(path_xml_file)
	root = tree.getroot()
	#
	metabolites=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}species')
	genes=root.findall('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct')
	reactions=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}reaction')
	#
	#subsystems=[p.text.split(': ')[1] for p in root.findall('.//{http://www.w3.org/1999/xhtml}p') if 'SUBSYSTEM' in p.text]
	#subsystems=list(set([subsystem.split('|')[1] for subsystem in subsystems]))
	#
	details=[len(genes),len(reactions),len(metabolites)]
	model_stats.loc[oid]=details
	#
	#
	#
	# get all reactions_to_subsystem dictionary
	dict_subsystems={reaction:p.text.split('|')[1] for reaction in reactions for p in reaction.findall('.//{http://www.w3.org/1999/xhtml}p') if 'SUBSYSTEM' in p.text}
	details_rxns=[]
	details_genes=[]
	for subsystem in subsystems:
	#for subsystem in model_stats_subsystem.columns:
		reactions_by_subsystem=[k for k,v in dict_subsystems.items() if v==subsystem]
		details_rxns.append(len(reactions_by_subsystem))
		genes_by_subsystem=set([gene.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct'] for reaction in reactions_by_subsystem for gene in reaction.findall('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProductRef')])
		details_genes.append(len(genes_by_subsystem))
	#
	model_stats_subsystem_rxns.loc[oid]=details_rxns
	model_stats_subsystem_genes.loc[oid]=details_genes



model_stats_subsystem_genes.to_csv('lookover.txt',sep='\t')
open('subsystems_to_curate.txt','w').write('\n'.join(subsystems))


beta=[v for k,v in dict_subsystems.items() if v==subsystem]


model_stats_subsystem_genes.to_csv('subsystems.txt',sep='\t')

organisms=organisms.drop([6,21])

################################################

page=open('evolview_model_stats_template.txt','r').read()
#open('evolview_model_stats_genes')

page=page.replace('Brettanomyces','Dekkera').replace('Hypocrea jecorina','Trichoderma reesei')

for order,row in organisms.iterrows():
	oid=organisms['oid'][order]
	print oid
	name=oid+' | '+organisms['Species'][order]
	line=name+'\t'+','.join([str(x) for x in model_stats.loc[oid].tolist()[:3]])
	open('evolview_model_stats.txt','a').write(line+'\n')



################################################

df_subsystem=pd.read_csv('subsystem_classification.txt',sep='\t',header=None).fillna('')

organisms=pd.read_excel('/Volumes/5TB/Organized/github/aybrah_old/aybrah_0.1.12_production.xlsx',sheet_name='curated_taxonomy_fungi').fillna('').set_index('order')
organisms=organisms.drop([6,21])


model_stats_subsystem_genes['other_total']=''
model_stats_subsystem_rxns['other_total']=''

for order,row in organisms.iterrows():
	oid=organisms['oid'][order]
	print oid
	subsystems_other=[df_subsystem[0][index] for index,row in df_subsystem[df_subsystem[1]=='other'].iterrows()]
	#
	rxns_other=sum([model_stats_subsystem_rxns[subsystem_other][oid] for subsystem_other in subsystems_other])
	genes_other=sum([model_stats_subsystem_genes[subsystem_other][oid] for subsystem_other in subsystems_other])
	model_stats_subsystem_rxns['other_total'][oid]=rxns_other
	model_stats_subsystem_genes['other_total'][oid]=genes_other


colors=['#a6cee3',
'#1f78b4',
'#b2df8a',
'#33a02c',
'#fb9a99',
'#e31a1c',
'#fdbf6f',
'#ff7f00',
'#cab2d6',
'#6a3d9a',
'#ffff99']


colors=['#8dd3c7',
#'#ffffb3',
'#bebada',
'#fb8072',
'#80b1d3',
'#fdb462',
'#b3de69',
'#fccde5',
'#999999',
'#bc80bd']
#'#ccebc5']
#'#ffed6f']


for order,row in organisms.iterrows():
	oid=organisms['oid'][order]
	print oid
	name=oid+' | '+organisms['Species'][order]
	#
	subsystems_nonother=[df_subsystem[0][index] for index,row in df_subsystem[df_subsystem[1]!='other'].iterrows()]
	subsystems_process=subsystems_nonother+['other_total']
	rxns_string=','.join([str(model_stats_subsystem_rxns[subsystem][oid]) for subsystem in subsystems_process])
	genes_string=','.join([str(model_stats_subsystem_genes[subsystem][oid]) for subsystem in subsystems_process])
	line_rxns=name+'\t'+rxns_string
	line_genes=name+'\t'+genes_string
	open('evolview_model_stats_rxns.txt','a').write(line_rxns+'\n')
	#open('evolview_model_stats_genes.txt','a').write(line_genes+'\n')


output_colors=[]

for i,subsystem in enumerate(subsystems_process):
	i,subsystem
	output_colors.append(colors[i%len(colors)])

page=open('evolview_model_stats_genes.txt','r').read()
page=page.replace('Brettanomyces','Dekkera').replace('Hypocrea jecorina','Trichoderma reesei')
page=page.replace('$colors$',','.join(output_colors))
page=page.replace('$title$','Titles')
page=page.replace('$dots$',','.join(subsystems_process))
open('evolview_model_stats_genes4.txt','w').write(page)

########################################################################


output_colors=[]

for i,subsystem in enumerate(subsystems_process):
	i,subsystem
	output_colors.append(colors[i%len(colors)])

page=open('evolview_model_stats_rxns.txt','r').read()
page=page.replace('Brettanomyces','Dekkera').replace('Hypocrea jecorina','Trichoderma reesei')
page=page.replace('$colors$',','.join(output_colors))
page=page.replace('$title$','Titles')
page=page.replace('$dots$',','.join(subsystems_process))
open('evolview_model_stats_rxns4.txt','w').write(page)


yields=pd.read_excel('strain_demand_yields_gram_big.xlsx')
for order,row in organisms.iterrows():
	oid=organisms['oid'][order]
	name=oid+' | '+organisms['Species'][order]
	output=','.join([str(yields[d][oid]) for d in debugs[5:25]])
	line=name+'\t'+output
	open('heat_map.txt','a').write(line+'\n')

#rgm | Rhodotorula graminis	text=Pucciniomycotina,fontcolor=black,linewidth=4,bkcolor=#3288bd


page=open('tree_fungal.newick','r').read()

for index,row in organisms.iterrows():
	oid=organisms.oid[index]
	name=organisms.Species[index]
	print name
	readable=taxonomy[(taxonomy.oids==oid) & (taxonomy.level=='Strain')]['readable'].item()
	output=oid+' | '+name+'\ttext='+readable+',fontcolor=black,linewidth=4,bkcolor=white'
	open('evolview_model_name.txt','a').write(output+'\n')



