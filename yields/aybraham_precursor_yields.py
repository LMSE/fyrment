

"""

test the yields of metabolic network

"""

import pandas as pd
import cobra
import urllib

url_aybrah_xlsx='https://github.com/kcorreia/aybrah/raw/master/aybrah.xlsx'

taxonomy=pd.read_excel(url_aybrah_xlsx,sheet_name='taxon_nodes')
#taxonomy.drop(46)
drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)


taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))[1]




version_aybrah=urllib.urlopen('https://github.com/kcorreia/aybrah/raw/master/version.txt').read()
version_aybraham=urllib.urlopen('https://github.com/kcorreia/aybraham/raw/master/version.txt').read()




organisms=pd.read_excel(url_aybrah_xlsx,sheet_name='curated_taxonomy_fungi').fillna('').set_index('order')
organisms=organisms.drop([6,21])


debugs=urllib.urlopen('https://github.com/kcorreia/aybraham/raw/master/version.txt').read().split('\n')
debugs=filter(None,debugs)
#debugs=filter(None,open('debug_rxns.txt','r').read().split('\n'))


oid_to_rxns={}

df_temp=pd.DataFrame(columns=debugs+['ATPM'])




#for order,row in organisms.iterrows():
for index,row in taxonomy[taxonomy['level']=='Strain'].iterrows():
	name=taxonomy['name'][index]
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	readable=taxonomy['readable'][index]
	path_xml_dir = '../genre/xml/'+str(level_order).zfill(2)+'_'+level.lower()
	path_xml_file=path_xml_dir+'/i'+readable+'.xml'
	oid=taxonomy['oids'][index]
	if level=='Strain':
		print name
	#
	model=cobra.io.read_sbml_model(path_xml_file)
	# get the maximum values
	output=[]
	for debug in debugs:
		model.reactions.get_by_id(debug).bounds=(0.0,1000.0)
		model.objective=debug
		fba=model.optimize()
		value=round(model.objective.value,4)
		output.append(value)
		fluxes_active=[rxnid for rxnid in fba.fluxes.keys() if abs(fba.fluxes[rxnid])>0.01]
		open('aybraham_'+version_aybraham+'_fluxes_active.txt','a').write('\t'.join([oid,debug,';'.join(fluxes_active)])+'\n')
		model.reactions.get_by_id(debug).bounds=(0.0,0.0)
	#
	model.objective='ATPM'
	fba=model.optimize()
	'ATPM'
	value=round(model.objective.value,4)
	output.append(value)
	fluxes_active=[rxnid for rxnid in fba.fluxes.keys() if abs(fba.fluxes[rxnid])>0.01]
	open('aybraham_'+version_aybraham+'_fluxes_active.txt','a').write('\t'.join([oid,'ATPM',';'.join(fluxes_active)])+'\n')
	#
	df_temp.loc[oid]=output
	oid_to_rxns[oid]=[reaction.id for reaction in model.reactions]

#cobra.io.sbml3.validate_sbml_model('./xml/'+model_xml)

df_temp.to_csv('aybraham_'+version_aybraham+'_yields.txt',sep='\t')




############################################################


"""
tolerance=0.01
yields=pd.read_csv('aybraham_'+version_aybraham+'_yields.txt',sep='\t').set_index('Unnamed: 0')
#organisms=pd.read_excel('/Volumes/5TB/Organized/github/aybrah_old/aybrah_'+version_aybrah+'_production.xlsx',sheet_name='curated_taxonomy_fungi').fillna('')


fluxes_active=pd.read_csv('aybraham_'+version_aybraham+'_fluxes_active.txt',sep='\t',header=None)

fluxes_active['sai']=''
fluxes_active['pic']=''
fluxes_active['yli']=''
fluxes_active['zro']=''

for oid,row in yields.iterrows():
	oid
	gapfilling_needed=[objective for objective in yields.columns if abs(yields[objective][oid])<tolerance]
	#gapfilling_needed.remove('DEBUG_PROTEIN')
	for debug in gapfilling_needed:
		debug
		rxns_debug_pic=fluxes_active[(fluxes_active[1]==debug) & (fluxes_active[0]=='pic')][2].item().split(';')
		rxns_debug_zro=fluxes_active[(fluxes_active[1]==debug) & (fluxes_active[0]=='zro')][2].item().split(';')
		rxns_debug_sai=fluxes_active[(fluxes_active[1]==debug) & (fluxes_active[0]=='sai')][2].item().split(';')
		rxns_debug_yli=fluxes_active[(fluxes_active[1]==debug) & (fluxes_active[0]=='yli')][2].item().split(';')
		#
		alpha=set(rxns_debug_sai)-set(oid_to_rxns[oid])
		beta=set(rxns_debug_pic)-set(oid_to_rxns[oid])
		gamma=set(rxns_debug_yli)-set(oid_to_rxns[oid])
		delta=set(rxns_debug_zro)-set(oid_to_rxns[oid])
		index=fluxes_active[(fluxes_active[0]==oid) & (fluxes_active[1]==debug)].index.item()
		fluxes_active['sai'][index]=';'.join(alpha)
		fluxes_active['pic'][index]=';'.join(beta)
		fluxes_active['yli'][index]=';'.join(gamma)
		fluxes_active['zro'][index]=';'.join(delta)



fluxes_active.to_csv('aybraham_'+version_aybraham+'_fluxes_active.txt',sep='\t',index=False)

"""

