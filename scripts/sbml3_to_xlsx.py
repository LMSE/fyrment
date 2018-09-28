

import pandas as pd
import xml.etree.ElementTree as ET
from pandas import ExcelWriter
import os



def write_xml_to_xlsx(level_order,level,readable):
	path_xml = './genre/xml/'+str(level_order).zfill(2)+'_'+level.lower()+'/i'+readable+'.xml'
	tree = ET.parse(path_xml)
	root = tree.getroot()
	#
	dict_parameters={parameter.attrib['id']:parameter.attrib['value'] for parameter in root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}parameter')}
	#
	genes=root.findall('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct')
	metabolites=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}species')
	compartments=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}compartment')
	reactions=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}reaction')
	metabolites=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}species')
	objectives=root.findall('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}listOfObjectives')
	#
	objective_dict={reaction_objective.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}reaction'][2:]: \
	reaction_objective.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}coefficient']\
		for objective in objectives \
		for reaction_objective in objective.findall('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}fluxObjective')}
	#
	names='Rxn name|Rxn description|Formula|Gene-reaction association|Genes|Protein|Subsystem|Reversible|LB|UB|Objective|Notes|References'.split('|')
	model_rxns=pd.DataFrame(columns=names)
	#
	for reaction in reactions:
		rxn_name=reaction.attrib['id'][2:]
		#if rxn_name=='HEX1':
		#	break
		print rxn_name
		rxn_description=reaction.attrib['name']
		ub=dict_parameters[reaction.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}upperFluxBound']]
		lb=dict_parameters[reaction.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound']]
		subsystem=[p.text.split(': ')[1] for p in reaction.findall('.//{http://www.w3.org/1999/xhtml}p') if 'SUBSYSTEM' in p.text][0]
		# reactants
		print('\tparse formula')
		reactants=[species for list_of_reactants in reaction.findall('.//{http://www.sbml.org/sbml/level3/version1/core}listOfReactants')\
			for species in list_of_reactants.findall('.//{http://www.sbml.org/sbml/level3/version1/core}speciesReference')]
		reactants_readable=[species.attrib['stoichiometry']+' '+species.attrib['species'][2:-2]+'['+species.attrib['species'][-1:]+']' \
					if species.attrib['stoichiometry']!='1'\
					 else species.attrib['species'][2:-2]+'['+species.attrib['species'][-1:]+']' for species in reactants]
		# products
		products=[species for list_of_reactants in reaction.findall('.//{http://www.sbml.org/sbml/level3/version1/core}listOfProducts')\
			for species in list_of_reactants.findall('.//{http://www.sbml.org/sbml/level3/version1/core}speciesReference')]
		products_readable=[species.attrib['stoichiometry']+' '+species.attrib['species'][2:-2]+'['+species.attrib['species'][-1:]+']' \
					if species.attrib['stoichiometry']!='1'\
					 else species.attrib['species'][2:-2]+'['+species.attrib['species'][-1:]+']' for species in products]
		print('\treversibility')
		if reaction.attrib['reversible']=='true':
			arrow=' <=> '
			reversible=1
		else:
			arrow=' -> '
			reversible=0
		if rxn_name in objective_dict.keys():
			objective_value=1
		else:
			objective_value=0
		formula=' + '.join(reactants_readable)+arrow+' + '.join(products_readable)
		# multiple complexes, isozymes
		gprs=reaction.find('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProductAssociation')
		output=[]
		complex_type=''
		if gprs is not None:
			for gpr in gprs.getchildren():
				gpr
				# isozymes??
				if gpr.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProductRef':
					string=gpr.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct']
					#string=string[2:]
					output.append(string)
				# isozymes (could be single enzymes or complexes)
				if gpr.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}or':
					for variant in gpr.getchildren():
						if variant.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}and':
							# includes nested paralogs with the OR operator
							string=[gene.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct'] if gene.tag!='{http://www.sbml.org/sbml/level3/version1/fbc/version2}or'\
															else '( '+'  |  '.join([ g.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct']  for g in gene.getchildren()  ])+' )' for gene in variant.getchildren()]
							#if level=='Strain':
								#string=[s[2:] for s in string]
							output.append('( '+'  &  '.join(string)+' )')
						# isozyme??
						elif variant.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProductRef':
							string=variant.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct']
							#if level=='Strain':
							#	string=string[2:]
							output.append(string)
				# complex, no multiple variants
				if gpr.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}and':
					output_complex=[]
					for variant in gpr.getchildren():
						if variant.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}or':
							string=[gene.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct'] for gene in variant.findall('.//{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProductRef')]
							#if level=='Strain':
							#	string=[s[2:] for s in string]
							output_complex.append('( '+'  |  '.join(string)+' )')
						elif variant.tag=='{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProductRef':
							string=variant.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}geneProduct']
							#if level=='Strain':
							#	string=string[2:]
							output_complex.append(string)
					output.append('  &  '.join(output_complex))
		# should specify || or | ?
		gpr=' | '.join(output)
		notes=''.join([p.text.split('NOTES: ')[1] for p in reaction.find('.//{http://www.sbml.org/sbml/level3/version1/core}notes').findall('.//{http://www.w3.org/1999/xhtml}p') if 'NOTES: ' in p.text])
		# need to populate here
		references=''
		details=[rxn_name,rxn_description,formula,gpr,'','',subsystem,reversible,lb,ub,objective_value,notes,references]
		model_rxns.loc[len(model_rxns)]=details
	#
	names='Metabolite name|Metabolite description|Metabolite neutral formula|Metabolite charged formula|Metabolite charge|Metabolite Compartment|Metabolite KEGG ID|Metabolite PubChem ID|Metabolite CheBI ID|InChI|SMILES'.split('|')
	model_mets=pd.DataFrame(columns=names)
	#
	for metabolite in metabolites:
		#
		met_id=metabolite.attrib['id'][2:-2]+'['+metabolite.attrib['id'][-1:]+']'
		print met_id
		met_description=metabolite.attrib['name']
		if '{http://www.sbml.org/sbml/level3/version1/fbc/version2}charge' in metabolite.attrib.keys():
			charge=metabolite.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}charge']
		else:
			charge=''
		if '{http://www.sbml.org/sbml/level3/version1/fbc/version2}chemicalFormula' in metabolite.attrib.keys():
			formula_charged=metabolite.attrib['{http://www.sbml.org/sbml/level3/version1/fbc/version2}chemicalFormula']
		else:
			formula_charged=''
		compartment=metabolite.attrib['id'].split('_')[-1]
		#
		output=[reference.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'].split('identifiers.org/')[1] for reference in metabolite.findall('.//{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li')]
		dict_metabolite={o.split('/')[0]:'/'.join(o.split('/')[1:]) for o in output}
		#
		if 'inchi' in dict_metabolite.keys():
			inchi=dict_metabolite['inchi']
		else:
			inchi=''
		#
		if 'pubchem.compound' in dict_metabolite.keys():
			pubchem_compound=dict_metabolite['pubchem.compound']
		else:
			pubchem_compound=''
		#
		if 'chebi' in dict_metabolite.keys():
			chebi=dict_metabolite['chebi']
		else:
			chebi=''
		#
		if 'kegg.compound' in dict_metabolite.keys():
			kegg=dict_metabolite['kegg.compound']
		else:
			kegg=''
		#
		if 'metacyc.compound' in dict_metabolite.keys():
			metacyc=dict_metabolite['metacyc.compound']
		else:
			metacyc=''
		smiles=''
		details=[met_id,met_description,'',formula,charge,compartment,kegg,pubchem_compound,chebi,inchi,smiles]
		model_mets.loc[len(model_mets)]=details
	path_xlsx_dir = './genre/xlsx/'+str(level_order).zfill(2)+'_'+level.lower()+'/'
	if not os.path.exists(path_xlsx_dir):
		os.makedirs(path_xlsx_dir)
	path_xlsx_file = './genre/xlsx/'+str(level_order).zfill(2)+'_'+level.lower()+'/i'+readable+'.xlsx'
	writer = ExcelWriter(path_xlsx_file)
	model_rxns.to_excel(writer,'reactions',index=False)
	model_mets.to_excel(writer,'metabolites',index=False)
	writer.save()




url_aybrah_xlsx='https://github.com/kcorreia/aybrah/raw/master/aybrah.xlsx'

taxonomy=pd.read_excel(url_aybrah_xlsx,sheet_name='taxon_nodes')
#taxonomy.drop(46)
drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)


taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))[1]


for index,row in taxonomy.iterrows():
#for index,row in taxonomy[taxonomy.level=='Strain'].iterrows():
#for index in [111]:
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	name=taxonomy['name'][index]
	print name
	path = './genre/xlsx/'+str(level_order).zfill(2)+'_'+level.lower()
	readable=taxonomy['readable'][index].replace('.','_')
	write_xml_to_xlsx(level_order,level,readable)

"""

model_xml='13_Strain_Scheffersomyces_stipitis_CBS_6054_AYbRAHAM_0_0_6_GENRE.xml'
tree = ET.parse('./xml/'+model_xml)
root = tree.getroot()


reactions=root.findall('.//{http://www.sbml.org/sbml/level3/version1/core}reaction')
for reaction in reactions:
	if reaction.attrib['id']=='R_GAPD':
		break

"""

