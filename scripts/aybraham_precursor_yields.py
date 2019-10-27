

"""

test the yields of metabolic network

"""

import pandas as pd
import cobra
import urllib

#url_aybrah_xlsx='https://github.com/kcorreia/aybrah/raw/master/aybrah.xlsx'
path_aybrah_xlsx='../aybrah/aybrah.xlsx'
taxonomy=pd.read_excel(path_aybrah_xlsx,sheet_name='taxon_nodes')
#taxonomy.drop(46)
drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)


taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))[1]




version_aybrah=urllib.urlopen('../aybrah/version.txt').read()

#version_aybraham=urllib.urlopen('https://github.com/kcorreia/aybraham/raw/master/version.txt').read()
version_fyrment=open('version.txt','r').read()



organisms=pd.read_excel(path_aybrah_xlsx,sheet_name='curated_taxonomy_fungi').fillna('').set_index('order')

indices_to_drop=organisms[organisms.oid.str.contains('ani|cpr')].index.tolist()
organisms=organisms.drop(indices_to_drop)


#debugs=urllib.urlopen('https://github.com/kcorreia/aybraham/raw/master/version.txt').read().split('\n')
#debugs=filter(None,debugs)
debugs=filter(None,open('./yields/debug_rxns.txt','r').read().split('\n'))


oid_to_rxns={}

df_temp=pd.DataFrame(columns=debugs+['ATPM'])

"""
109
125

rxns_pku=[r.id for r in model.reactions]

rxns_pme=[r.id for r in model.reactions]
for x in list(set(rxns_pme)-set(rxns_pku)):
	print(x)
"""

#for order,row in organisms.iterrows():
for index,row in taxonomy[taxonomy['level']=='Strain'].iterrows():
	name=taxonomy['name'][index]
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	readable=taxonomy['readable'][index]
	path_xml_dir = './genre/xml/'+str(level_order).zfill(2)+'_'+level.lower()
	path_xml_file=path_xml_dir+'/i'+readable+'.xml'
	oid=taxonomy['oids'][index]
	#cobra.io.sbml3.validate_sbml_model(path_xml_file)
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
		open('aybraham_'+version_fyrment+'_fluxes_active.txt','a').write('\t'.join([oid,debug,';'.join(fluxes_active)])+'\n')
		model.reactions.get_by_id(debug).bounds=(0.0,0.0)
	#
	model.objective='ATPM'
	fba=model.optimize()
	'ATPM'
	value=round(model.objective.value,4)
	output.append(value)
	fluxes_active=[rxnid for rxnid in fba.fluxes.keys() if abs(fba.fluxes[rxnid])>0.01]
	open('aybraham_'+version_fyrment+'_fluxes_active.txt','a').write('\t'.join([oid,'ATPM',';'.join(fluxes_active)])+'\n')
	#
	df_temp.loc[oid]=output
	oid_to_rxns[oid]=[reaction.id for reaction in model.reactions]



#cobra.io.sbml3.validate_sbml_model('./xml/'+model_xml)

df_temp.to_csv('./yields/aybraham_'+version_fyrment+'_yields.txt',sep='\t')

df_temp=pd.read_csv('./yields/aybraham_'+version_fyrment+'_yields.txt',sep='\t').set_index('Unnamed: 0')

def get_MW(aa):
	metabolites=model.reactions.get_by_id(aa).metabolites
	for metabolite in metabolites.keys():
		#print(aa,metabolites)
		return metabolite.formula_weight

for oid,row in df_temp.iterrows():
	#break
	index=taxonomy[(taxonomy.level=='Strain') & (taxonomy.oids==oid)].index.item()
	readable=taxonomy['readable'][index]
	name=' '.join(taxonomy.name[index].split(' ')[:2])
	row=oid+' | '+name+'\t'+','.join([str(round(df_temp[aa][oid]/10*get_MW(aa)/180,4)) for aa in df_temp.columns[5:-24]])
	print(row)







import cobra
from os import listdir
from os.path import isfile, join



import pandas as pd





def add_reactions():
	def check_met(metid_compartment):
		if metid_compartment not in model.metabolites:
			M_x = cobra.Metabolite(
		    	metid_compartment,
		    	formula='',
		    	name=metid_compartment,
		    	compartment=metid_compartment.split('_')[-1])
		else:
			M_x=model.metabolites.get_by_id(metid_compartment)
		return M_x
	#########
	"""
	Add reactions for metabolite maximization
	"""
	############################################
	reaction_bdo1 = cobra.Reaction('BDO1')
	reaction_bdo1.name = 'Aldolase'
	reaction_bdo1.subsystem = 'Aldolase-BDO pathway'
	reaction_bdo1.lower_bound = 0.  # This is the default
	reaction_bdo1.upper_bound = 1000.  # This is the default
	#
	M_3hbtal_c = cobra.Metabolite(
	    '3hbtal_c',
	    formula='',
	    name='3-hydroxybutanal',
	    compartment='c')
	#
	M_acald_c=model.metabolites.get_by_id('acald_c')
	#
	reaction_bdo1.add_metabolites({
	    M_acald_c: -2.0,
	    M_3hbtal_c: 1.0
	})
	reaction_bdo2 = cobra.Reaction('BDO2')
	reaction_bdo2.name = 'Aldolase-BDO pathway'
	reaction_bdo2.subsystem = '13BDO pathway'
	reaction_bdo2.lower_bound = 0.  # This is the default
	reaction_bdo2.upper_bound = 1000.  # This is the default
	#
	M_13bdo__R_c = cobra.Metabolite(
	    '13bdo__R_c',
	    formula='',
	    name='1,3-butandiol',
	    compartment='c')
	#
	M_nadph_c=model.metabolites.get_by_id('nadph_c')
	M_nadp_c=model.metabolites.get_by_id('nadp_c')
	M_h_c=model.metabolites.get_by_id('h_c')
	#
	reaction_bdo2.add_metabolites({
	    M_3hbtal_c: -1.0,
	    M_nadph_c: -1.0,
	    M_h_c: -1.0,
	    M_nadp_c: 1.0,
	    M_13bdo__R_c: 1.0
	})
	############################################
	reaction_bdo3 = cobra.Reaction('DM_13bdo__R')
	reaction_bdo3.name = 'Aldolase-BDO pathway'
	reaction_bdo3.subsystem = '13BDO pathway'
	reaction_bdo3.lower_bound = 0.  # This is the default
	reaction_bdo3.upper_bound = 1000.  # This is the default
	#
	reaction_bdo3.add_metabolites({
	    M_13bdo__R_c: -1.0
	})
	#
	model.add_reactions([reaction_bdo1, reaction_bdo2, reaction_bdo3])
	############################################
	M_lycop_c=check_met('lycop_c')
	#if 'lycop_c' in model.metabolites:
	#	M_lycop_c=model.metabolites.get_by_id('lycop_c')
	#else:
	#	M_lycop_c = cobra.Metabolite(
	#   	'lycop_c',
	#    	formula='',
	#    	name='lycopene',
	#    	compartment='c')
	#
	reaction_EX_lycop = cobra.Reaction('DM_lycop_c')
	reaction_EX_lycop.name = 'Demand lycopene'
	reaction_EX_lycop.subsystem = 'Exchange'
	reaction_EX_lycop.lower_bound = 0.  # This is the default
	reaction_EX_lycop.upper_bound = 1000.  # This is the default
	#
	reaction_EX_lycop.add_metabolites({
	    M_lycop_c: -1.0
	})
	#
	model.add_reactions([reaction_EX_lycop])
	############################################
	M_ptrc_c=check_met('ptrc_c')
	reaction_EX_ptrc = cobra.Reaction('DM_ptrc_c')
	reaction_EX_ptrc.name = 'Demand ptrc'
	reaction_EX_ptrc.subsystem = 'Demand'
	reaction_EX_ptrc.lower_bound = 0.  # This is the default
	reaction_EX_ptrc.upper_bound = 1000.  # This is the default
	#
	reaction_EX_ptrc.add_metabolites({
	    M_ptrc_c: -1.0
	})
	model.add_reactions([reaction_EX_ptrc])
	############################################
	M_msa_m=check_met('msa_m')
	#if 'msa_m' not in model.metabolites:
	#	M_msa_m = cobra.Metabolite(
	#    	'msa_m',
	#    	formula='',
	#    	name='msa_m',
	#    	compartment='c')
	#else:
	#	M_msa_m=model.metabolites.get_by_id('msa_m')
	M_msa_c=model.metabolites.get_by_id('msa_c')
	#
	reaction_TR_msa = cobra.Reaction('TR_msa')
	reaction_TR_msa.name = 'msa transport'
	reaction_TR_msa.subsystem = 'Transport'
	reaction_TR_msa.lower_bound = -1000.  # This is the default
	reaction_TR_msa.upper_bound = 1000.  # This is the default
	#
	reaction_TR_msa.add_metabolites({
	    M_msa_m: -1.0,
	    M_msa_c: 1.0
	})
	model.add_reactions([reaction_TR_msa])
	############################################
	M_frdp_c=model.metabolites.get_by_id('frdp_c')
	M_ppi_c=model.metabolites.get_by_id('ppi_c')
	M_farnesene_c=check_met('farnesene_c')
	reaction_farnesene_syn = cobra.Reaction('farnesene_syn')
	reaction_farnesene_syn.name = 'farnesene synthase'
	reaction_farnesene_syn.subsystem = 'farnesene synthase'
	reaction_farnesene_syn.lower_bound = 0.  # This is the default
	reaction_farnesene_syn.upper_bound = 1000.  # This is the default
	#
	reaction_farnesene_syn.add_metabolites({
	    M_frdp_c: -1.0,
	    M_farnesene_c: 1.0,
	    M_ppi_c:	1.0
	})
	reaction_DM_farnesene = cobra.Reaction('DM_farnesene')
	reaction_DM_farnesene.name = 'Farnesene demand reaction'
	reaction_DM_farnesene.subsystem = 'Demand reaction'
	reaction_DM_farnesene.lower_bound = 0.  # This is the default
	reaction_DM_farnesene.upper_bound = 1000.  # This is the default
	#
	reaction_DM_farnesene.add_metabolites({
	    M_farnesene_c: -1.0,
	})
	model.add_reactions([reaction_farnesene_syn, reaction_DM_farnesene])
	############################################
	M_3dhsk_c=model.metabolites.get_by_id('3dhsk_c')
	M_catechol_c=model.metabolites.get_by_id('catechol_c')
	M_co2_c=model.metabolites.get_by_id('co2_c')
	M_nadh_c=model.metabolites.get_by_id('nadh_c')
	M_nad_c=model.metabolites.get_by_id('nad_c')
	M_h_c=model.metabolites.get_by_id('h_c')
	M_h2o_c=model.metabolites.get_by_id('h2o_c')
	M_adipate_c=check_met('adipate_c')
	M_ccmuac_c=check_met('ccmuac_c')
	M_protocatechuate_c=check_met('protocatechuate_c')
	#
	adipic_acid1 = cobra.Reaction('DHSD_adipate')
	adipic_acid1.name = 'DHS-D'
	adipic_acid1.subsystem = 'adipic_acid'
	adipic_acid1.lower_bound = 0.  # This is the default
	adipic_acid1.upper_bound = 1000.  # This is the default
	adipic_acid1.add_metabolites({
	    M_3dhsk_c: -1.0,
	    M_protocatechuate_c: 1.0,
	    M_h2o_c:	1.0
	})
	#
	adipic_acid2 = cobra.Reaction('PCAD_adipate')
	adipic_acid2.name = 'PCA-D'
	adipic_acid2.subsystem = 'adipic_acid'
	adipic_acid2.lower_bound = 0.  # This is the default
	adipic_acid2.upper_bound = 1000.  # This is the default
	adipic_acid2.add_metabolites({
	    M_protocatechuate_c: -1.0,
	    M_h_c:-1.0,
	    M_catechol_c: 1.0,
	    M_co2_c:	1.0
	})
	#
	adipic_acid3 = cobra.Reaction('ER_adipate')
	adipic_acid3.name = 'ER'
	adipic_acid3.subsystem = 'adipic_acid'
	adipic_acid3.lower_bound = 0.  # This is the default
	adipic_acid3.upper_bound = 1000.  # This is the default
	adipic_acid3.add_metabolites({
	    M_ccmuac_c: -1.0,
	    M_nadh_c: -2.0,
	    M_h_c:		-2.0,
	    M_nad_c:	2.0,
	    M_adipate_c:	1.0
	})
	adipic_acid4 = cobra.Reaction('DM_adipate')
	adipic_acid4.name = 'DM'
	adipic_acid4.subsystem = 'adipic_acid'
	adipic_acid4.lower_bound = 0.  # This is the default
	adipic_acid4.upper_bound = 1000.  # This is the default
	adipic_acid4.add_metabolites({
	    M_adipate_c: -1.0
	})
	model.add_reactions([adipic_acid1,adipic_acid2,adipic_acid3,adipic_acid4])
	############################################
	M_3hpcoa_x=check_met('3hpcoa_x')
	M_3hpcoa_m=check_met('3hpcoa_m')
	#
	msa_precursor = cobra.Reaction('3HPCOA_TR')
	msa_precursor.name = 'msa transport'
	msa_precursor.subsystem = 'msa'
	msa_precursor.lower_bound = 0.  # This is the default
	msa_precursor.upper_bound = 1000.  # This is the default
	msa_precursor.add_metabolites({
	    M_3hpcoa_x: -1.0,
	    M_3hpcoa_m:	1.0
	})





def block_rxns():
	for rxnid in 'ALCD24xi|ALCD22xi|ALCD25xi|TYRSLDHx|ALCD26xi'.split('|'):
		if rxnid in model.reactions:
			model.reactions.get_by_id(rxnid).bounds=(0,0)




taxons=['12_strain', '04_subdivision', '00_domain']



targets=['EX_succ_e',
'EX_glyclt_e',
#'EX_msa_e',
#'EX_catechol_e',
'DM_ptrc_c',
# alcohols
'EX_iamoh_e',
'EX_2mbtoh_e',
'EX_2phetoh_e',
'EX_tyrosol_e',
'EX_ind3eth_e',
'EX_ibutoh_e',
# isoprenoids?
'DM_ERGST',
'DM_lycop_c',
'DM_farnesene',
#'DM_adipate',
'DM_13bdo__R']


path='./genre/xml/12_strain/iOgapa_DL_1.xml.gz'
#model=cobra.io.read_sbml_model(path)


#cobra.io.sbml3.validate_sbml_model(path)


df=pd.DataFrame(columns=targets)
df2={}

for taxon in taxons:
	#path_xml_file='/Volumes/SD250GB/Organized/github/redo/aybraham/genre/xml/'+taxon+'/'
	path_xml_file='./genre/xml/'+taxon+'/'
	models = [f for f in listdir(path_xml_file) if isfile(join(path_xml_file, f)) if '.xml' in f]
	print('get model names')
	for fname_model in models:
		readable=fname_model.split('.')[0]
		df2[readable]=pd.DataFrame()
		print(readable)
		model=cobra.io.read_sbml_model(path_xml_file+fname_model)
		add_reactions()
		block_rxns()
		# metabolite targets
		output1=[]
		output2=[]
		for target in targets:
			print(target)
			if target=='DM_ERGST':
				model.reactions.get_by_id('DM_ERGST').bounds=(0,1000)
			model.objective=target
			# turn off NADH reaction
			# optimize
			fba_solution=model.optimize()
			#
			if target=='DM_ERGST':
				model.reactions.get_by_id('DM_ERGST').bounds=(0,0)
			output1.append(round(fba_solution.objective_value,2))
			output2.append(pd.Series(fba_solution.fluxes,name=target))
		df.loc[readable]=output1
		df2[readable]=pd.concat(output2,axis=1)
		df2[readable].to_csv('./yields/'+readable+'.txt',sep='\t')



df.to_csv('./yields/meteng_yields.txt',sep='\t')

"""
Print for EvolView
"""

mw=pd.read_csv('./yields/meteng_formula_weights.txt',sep=',',header=None).set_index(0)

df.index=[i[1:] for i in df.index.tolist()]

for index,row in taxonomy[taxonomy.level=='Strain'].iterrows():
	oid=taxonomy['oids'][index]
	name=' '.join(taxonomy['name'][index].split(' ')[:2])
	readable=taxonomy['readable'][index]
	if readable in df.index:
	#if 'i'+readable in df.index:
		row=[round(df[exchange][readable]/10*mw[1][exchange]/180,2) for exchange in df.columns]
		print(oid+' | '+name+'\t'+','.join([str(r) for r in row]))





for index,row in taxonomy[taxonomy.level=='Strain'].iterrows():
	#break
	oid=taxonomy['oids'][index]
	name=' '.join(taxonomy['name'][index].split(' ')[:2])
	readable=taxonomy['readable'][index]
	fname='./yields/i'+readable+'.txt'
	df=pd.read_csv(fname,sep='\t').set_index('Unnamed: 0')
	#if 'i'+readable in df.index:
	row=[round(df[exchange][exchange]/10*mw[1][exchange]/180,2) for exchange in df.columns]
	print(oid+' | '+name+'\t'+','.join([str(r) for r in row]))




phylos=['iEukaryota.txt',
'iTaphrino.txt',
'iPuccinio.txt',
'iPezizo.txt',
'iBudding.txt']

for phylo in phylos:
	#break
	fname='./yields/'+phylo
	df=pd.read_csv(fname,sep='\t').set_index('Unnamed: 0')
	#if 'i'+readable in df.index:
	row=[round(df[exchange][exchange]/10*mw[1][exchange]/180,2) for exchange in df.columns]
	print(phylo.split('.')[0]+'\t'+','.join([str(r) for r in row]))





'meteng_yields.txt'



"""
for col in df.columns:
	df[col]=df[col].astype(float)
"""


"""
track changes for maximum yielding model
"""

reference='iSacce_S288c'

for target in targets:
	tracked=pd.DataFrame(columns='rxns1|rxns2|model'.split('|'))
	#target='DM_13bdo__R'
	df_ref=pd.read_csv('./yields/'+reference+'.txt',sep='\t')
	rxns_ref=df_ref[abs(df_ref[target])>0.001]['Unnamed: 0'].tolist()
	#
	df_delta=df[df[target]>df[target][reference]]
	for model_id,row in df_delta.iterrows():
		#df[target]
		df_sub=pd.read_csv('./yields/'+model_id+'.txt',sep='\t')
		rxns_sub=df_sub[abs(df_sub[target])>0.001]['Unnamed: 0'].tolist()
		#
		difference_used=set(rxns_sub)-set(rxns_ref)
		difference_present=difference_used-set(df_ref['Unnamed: 0'])
		difference_latent=difference_used-difference_present
		place_rxn(difference_present,difference_latent,model_id)
	tracked.to_csv(target+'.txt',sep='\t',index=False)


sorted(set(rxns_sub))

for x in :
	if x not in set(rxns_ref):
		breakker


def place_rxn(difference_present,difference_latent,model_id):
	text2=';'.join(sorted(difference_present))
	text1=';'.join(sorted(difference_latent))
	if text2 not in tracked['rxns2']:
		tracked.loc[len(tracked)]=[text1,text2,model_id]
	else:
		index=tracked[tracked['rxns2']==text2].index.item()
		tracked['model'][index]+';'+model_id



cobra.io.sbml3.validate_sbml_model(path_xml_file+fname_model)

#model.reactions.get_by_id('ALCD23xi').bounds=(0,0)
#model.reactions.get_by_id('GAPDi_nadp').bounds=(0,0)
#model.reactions.get_by_id('NADPPPS').bounds=(0,1000)

rxns='EX_succ_e|EX_glyclt_e|DM_ptrc_c|EX_iamoh_e|EX_2mbtoh_e|EX_2phetoh_e|EX_tyrosol_e|EX_ind3eth_e|EX_ibutoh_e|DM_ERGST|DM_lycop_c|DM_farnesene|DM_13bdo__R'.split('|')


aybrahm_organisms=pd.read_excel('https://github.com/kcorreia/aybraham/raw/master/aybraham.xlsx')


level

molar=pd.read_csv('meteng_compounds.txt',sep='\t').set_index('Reaction')

levels=['Strain','Subdivision','Domain']

for level in levels:
	for index,row in taxonomy[taxonomy.level==level].iterrows():
		readable=taxonomy['readable'][index]
		oid=taxonomy['oids'][index]
		name=taxonomy['name'][index]
		#name=organisms[organisms.oid==oid]['Species'].item()
		if oid in ['cpr','ani']:
			continue
		#elif level == 'Strain':
		#	name=organisms[organisms.oid==oid]['Strain'].item()
		if level=='Strain':
			name=organisms[organisms.oid==oid]['Species'].item()
			key=oid+' | '+name
		else:
			key=name
		output=key+'\t'+','.join([str(    round(df[rxn]['i'+readable]*molar['Molar'][rxn]/10/180,2   )) for rxn in rxns])+'\n'
		open('meteng_heat_map.txt','a').write(output)








for x in 'BZDH|PHEGLXDC|BZALDHy'.split('|'):
	fba_solution.fluxes[x]



model.objective='EX_13bdo__R'
model.reactions.get_by_id('EX_o2_e').bounds=(0,0)
fba_solution=model.optimize()









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

