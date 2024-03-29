


from libsbml import *
import pandas as pd
import urllib

import xml.etree.cElementTree as ET




"""

This script creates metabolic networks for different fungal taxonomies
by using FYRMENT and AYbRAH resources.


Example of a reaction in FYRMENT:
	rxnid: ALDD2y
	rxnd description: Aldehyde dehydrogenase (acetaldehyde, NADP)
	formula: 'acald[c] + nadp[c] + h2o[c] -> nadph[c] + ac[c] + 2 h[c]'
	FOG annotation: 'FOG00350  |  FOG00362  |  FOG00361  |  FOG00360'


Example of orthologous group annotation in AYbRAH:
	FOG			pta	sce	dbx	dha	spo pic
	FOG00350	0	0	0	0	1	0
	FOG00362	0	1	0	0	0	0
	FOG00361	1	0	1	0	0	0
	FOG00360	0	0	0	1	0	0


To determine whether a reaction is included, a reaction must have at least one
function protein/complex in the GPR (gene-protein-reaction) association.
| denotes the or operator, & denotes and operator used for complexes.
Sometimes proteins have different complex variants, in which case a double or operator is used ||:
	A & B || A & C
Complex variants with non-essential subunits are denoted with ?. 
In the example below D is non-essential,
which means that a GPR is included in a network if A and B are present.
	A & B & D?

Note that some reactions in FYRMENT are automatically included in all metabolic networks like:
BIOMASS, Exchange reactions, 


"""



def process_fog_annotation(oid,fog_annotation):
	"""
	Could and should clean this up.
	Purpose is to take in a fog_annotation in FYRMENT
	(such as 'FOG01654  |  FOG01655') and oid (such as kla, pic)
	to determine which  genes/proteins to include for the GPR.
	Checks if all essential subunits are present for protein complexes.
	Should not have nested functions in future versions
	"""
	def process_complex_variants(oid,fog_annotation):
		"""
		Code to process complex variants, which only looks at
		essential subunits. Examples include Complex I
		"""
		complex_check={}
		complex_variants=fog_annotation.split('||')
		for i,complex_variant in enumerate(complex_variants):
			complex_check[i]=[]
			subcomplexes=complex_variant.split('&&')
			for subcomplex in subcomplexes:
				subcomplex_list=subcomplex.replace('(','').replace(')','').replace(' ','').split('|')
				for subcomplex_component in subcomplex_list:
					fogs=subcomplex_component.split('&')
					for fog in fogs:
						complex_check[i].append(fog)
		complex_include={}
		for i in sorted(complex_check.keys()):
			# list of essential FOG's
			complex_essential=[fog[:8] for fog in complex_check[i] if '?' not in fog]
			# list of essential FOG's in oid
			complex_esssential_oid=[aybrah[oid][fog] for fog in complex_essential if len(aybrah[oid][fog[:8]])>0 ]
			# bypass if there is no oid (pan-model)
			#if len(inclusion_oids)>1:
			#	complex_include[i]=complex_check[i]
			# check if organism has all the designated essential proteins
			if len(complex_esssential_oid)==len(complex_essential):
				#oid=inclusion_oids[0]
				fogs_processed=[fog[:8] for fog in complex_check[i] if len(aybrah[oid][fog[:8]])>0]
				if len(fogs_processed)>0:
					complex_include[i]=fogs_processed
				print(oid,'include complex '+str(i))
			else:
				# return empty
				print('oid','not all subunits for complex '+str(i))
				#complex_include=None
		return complex_include
	# fog_annotation requires parsing
	if '&' in fog_annotation:
		# complex that meet essentialy requir.
		complex_include=process_complex_variants(oid,fog_annotation)
		fogs=complex_include
		#if oid != None:
		#	fogs=sorted(set([fog[:8] for key in complex_include for fog in complex_include[key] if len(aybra[oid][fog[:8]>0])]))
	# simple or fog_annotation
	else:
		fogs=filter(None,fog_annotation.replace('(','').replace(')','').replace(' ','').split('|'))
		if oid==None:
			fogs={i:[fog] for i,fog in enumerate(fogs)}
		else:
			fogs={i:[fog] for i,fog in enumerate(fogs) if len(aybrah[oid][fog])>0}
	return fogs



def parse_stoich_met(stoich_metabolite):
	"""
	Parse text formula to get stoich and metabolite
	"""
	# remove white space at the end of metabolite
	while stoich_metabolite[-1] in ' ':
		stoich_metabolite=stoich_metabolite[:-1]
	if ' ' in stoich_metabolite:
		stoich,met_id=stoich_metabolite.split(' ')
	else:
		stoich='1.0'
		met_id=str(stoich_metabolite)
	# necessary to remove features not compatible with SBML
	return stoich,met_id




def parse_formula(rxn_formula):
	"""
	Takes in rxn_formula to creates a dictionary of reactants and products
	with stoich/metabolite IDs.
	Example:
	Input:
	'6.0565 mannan[c] -> mannan_1g[c] + 5.0565 h2o[c]'
	Output:
	{'reactants': [('6.0565', 'mannan[c]')], 'products': [('1.0', 'mannan_1g[c]'),\
	 ('5.0565', 'h2o[c]')]}

	"""
	parsed_rxn={}
	parsed_rxn['reactants']=[]
	print(rxn_formula)
	if '->' in rxn_formula:
		# case of exchange reactions
		if len(list(filter(None,rxn_formula.split('->'))))==1:
			delim=' ->'
		# non-exchange reactions
		else:
			delim=' -> '
	elif '<=>' in rxn_formula:
		if len(list(filter(None,rxn_formula.split('<=>'))))==1:
			delim=' <=>'
		else:
			delim=' <=> '
	print('parse products')
	reactants=rxn_formula.split(delim)[0].split(' + ')
	for reactant in reactants:
		print('\t'+reactant)
		stoich,met_id=parse_stoich_met(reactant)
		parsed_rxn['reactants'].append((str(stoich),str(met_id)))
	# checks if products in the fomrula, combine with above conditions?
	if len(rxn_formula.split(delim))>1:
		print('parse products')
		parsed_rxn['products']=[]
		products=rxn_formula.split(delim)[1].split(' + ')
		for product in list(filter(None,products)):
			print('\t'+product)
			stoich,met_id=parse_stoich_met(product)
			parsed_rxn['products'].append((str(stoich),str(met_id)))
	return parsed_rxn


def get_rxn_inclusion():
	"""
	"""
	inclusion_automatic=['BIOMASS','DEMAND','DIFFUSION','EQUILIBRIUM','ESSENTIAL','EXCHANGE','SPONTANEOUS','ORPHAN']
	oids=aybrah.columns[2:].tolist()
	oid_rxns={oid:{} for oid in oids}
	for rxnid,row in rxns.iterrows():
		print(rxnid)
		#
		inclusion_model=rxns['model'][rxnid]
		inclusion_organism=rxns['organism'][rxnid].split('|')
		fog_annotation=rxns['FOG'][rxnid]
		# process the pan-reactome
		for oid in oids:
			#
			fogs=process_fog_annotation(oid,fog_annotation)
			# genetic evidence?
			if len(fogs)>0:
				# add index, rxnid as key for the reaction
				oid_rxns[oid][rxnid]='FOGS'
				#oid_rxns[oid][(index,rxnid)]='FOGS:'+'|'.join(fogs)
			# inclusion list
			elif inclusion_model in inclusion_automatic:
				# this will include genes for the ribosome even though essential genes are not fully annotated
				fogs=process_fog_annotation(oid,fog_annotation)
				if len(fogs)>0:
					oid_rxns[oid][rxnid]='FOGS'
				else:
					oid_rxns[oid][rxnid]=inclusion_model
			# organism specific evidence?
			elif oid in inclusion_organism:
				oid_rxns[oid][rxnid]='evidence in '+oid
			elif 'ORPHAN_FOG' in inclusion_model:
				fogs_check=inclusion_model.split(':')[1].split('|')
				fogs=[fog for fog in fogs_check if len(aybrah[oid][fog])>0 ]
				if len(fogs)>0:
					oid_rxns[oid][rxnid]='ORPHAN_FOG:'+'|'.join(fogs)
			elif 'ORPHAN_RXN' in inclusion_model:
				rxnid_checks=inclusion_model.split(':')[1].split('|')
				prereq_present=[ rxnid_check for rxnid_check in rxnid_checks if len(process_fog_annotation(oid,rxns['FOG'][rxnid_check]))>0]
				if len(prereq_present)>0:
					oid_rxns[oid][rxnid]='ORPHAN_RXN:'+'|'.join(prereq_present)
			elif 'PREREQ_RXN' in inclusion_model:
				rxnid_checks=inclusion_model.split(':')[1].split('|')
				prereq_present=[ rxnid_check for rxnid_check in rxnid_checks if len(process_fog_annotation(oid,rxns['FOG'][rxnid_check]))>0]
				if len(prereq_present)>0:
					oid_rxns[oid][rxnid]='PREREQ_RXN:'+'|'.join(prereq_present)
	return oid_rxns











def create_units():
	# Create unit definitions
	flux = model.createUnitDefinition()
	flux.setId('mmol_per_gDW_per_hr')
	#
	###################################
	unit = flux.createUnit()
	unit.setKind(UNIT_KIND_MOLE)
	unit.setExponent(1)
	unit.setScale(-3)
	unit.setMultiplier(1)
	###################################
	unit = flux.createUnit()
	unit.setKind(UNIT_KIND_GRAM)
	unit.setExponent(-1)
	unit.setScale(1)
	unit.setMultiplier(1)
	###################################
	unit = flux.createUnit()
	unit.setKind(UNIT_KIND_SECOND)
	unit.setExponent(-1)
	unit.setScale(1)
	unit.setMultiplier(1.0 * 60 * 60)



def create_compartments():
	###################################
	# Create compartments, name is set to id by default
	#
	for c,row in df_c.iterrows():
		c=str(c)
		print(c)
		c1 = model.createCompartment()
		status=c1.setId(c)
		check_status(status,'compartment_id',c)
		# needed to add CV term
		status=c1.setMetaId(c)
		status=c1.setName(str(df_c['compartment_name'][c]))
		check_status(status,'compartment_name',str(df_c['compartment_name'][c]))
		status=c1.setConstant(False)
		check_status(status,'compartment_id',c)
		#c1.setSize(1)
		#c1.setSpatialDimensions(3)
		#c1.setUnits('litre')
		value=df_c['GO term'][c]
		cv = CVTerm();
		cv.setQualifierType(BIOLOGICAL_QUALIFIER);
		cv.setBiologicalQualifierType(BQB_IS);
		cv.addResource('http://identifiers.org/go/'+str(value));
		c1.addCVTerm(cv);
		value_bigg=df_c['bigg.compartment'][c]
		if len(value_bigg)>0:
			cv = CVTerm();
			cv.setQualifierType(BIOLOGICAL_QUALIFIER);
			cv.setBiologicalQualifierType(BQB_IS);
			cv.addResource('http://identifiers.org/bigg.compartment/'+str(value_bigg));
			c1.addCVTerm(cv);






def get_metabolites_in_reactions():
	# gets all metabolites in the pan-model reaction network
	output_reaction_metabolites=[]
	for rxnid,row in rxns.iterrows():
		rxn_formula=rxns['Formula'][rxnid]
		if rxnid in inclusion_rxns:
			parsed_rxn=parse_formula(rxn_formula)
			output_reaction_metabolites.extend([met  for key in parsed_rxn.keys() for stoich,met in parsed_rxn[key] ])
	return output_reaction_metabolites












def create_species(name):
	# creates metabolic species with meta data (CHEBI, METACYC, etc)
	metabolites_in_reactions=get_metabolites_in_reactions()
	# Appears to be metabolites in pan-Fungi, but not captured in lower level GENRE
	if name=='Fungi':
		metabolites_in_metabolites=[met for met,row in mets.iterrows()]
		mets_to_add=[(met,met+'_'+c,c) for met,row in mets.iterrows()  for c in df_c.index.tolist() if met+'['+c+']' in metabolites_in_reactions]
		mets_not_in_reactions=list(set(metabolites_in_metabolites)-set([met[:-3] for met in metabolites_in_reactions]))
		mets_to_add+=[(met,met+'_c','c') for met in mets_not_in_reactions]
	else:
		mets_to_add=[(met,met+'_'+c,c) for met,row in mets.iterrows()  for c in df_c.index.tolist() if met+'['+c+']' in metabolites_in_reactions]
	for met,met_compartment,compartment in mets_to_add:
		#break
		print(met)
		name=str(mets['Metabolite description'][met])
		formula=str(mets['Metabolite charged formula'][met])
		charge=str(mets['Metabolite charge'][met])
		kegg_id=str(mets['Metabolite KEGG ID'][met])
		pubchem_cid=str(mets['Metabolite PubChem ID'][met])
		chebi_id=str(mets['Metabolite CheBI ID'][met])
		metacyc_id=str(mets['MetaCyc'][met])
		inchi=str(mets['InChI'][met])
		smiles=str(mets['SMILES'][met])
		inchikey=str(mets['InChIKey'][met])
		mnx=str(mets['MetNetX'][met])
		bigg=str(mets['BiGG'][met])
		# do not do anything with spider for now
		spider=str(mets['ChemSpiderID'][met])
		sbo=str(mets['SBO'][met])
		# assign to xml tag
		s1 = model.createSpecies()
		status=s1.setId(str('M_'+met_compartment))
		check_status(status,'species_id',str('M_'+met_compartment))
		status=s1.setMetaId('M_'+str(met_compartment))
		check_status(status,'species_metaid',str('M_'+met_compartment))
		status=s1.setName(name)
		check_status(status,'species_name',name)
		s1.setInitialConcentration(0)
		s1.setCompartment(str(compartment))
		s1.setConstant(False)
		s1.setBoundaryCondition(False)
		s1.setHasOnlySubstanceUnits(False)
		s1.setSBOTerm(sbo)
		#
		splugin = s1.getPlugin("fbc")
		if splugin != None: 
			if len(charge)>0:
				splugin.setCharge(int(float(charge)))
			if len(formula)>0:
				# chemical formula should follow the hill system
				splugin.setChemicalFormula(formula)
		#
		compound_identifiers=[
			('chebi',chebi_id),
			('kegg.compound',kegg_id),
			('metacyc.compound',metacyc_id),
			('pubchem.compound',pubchem_cid),
			('inchi',inchi),
			('inchikey',inchikey),
			('mnx.species',mnx),
			('bigg.species',bigg)]
		#
		for key,value in sorted(compound_identifiers):
			if len(value)>0:
				if key in ['pubchem.compound','chebi']:
					value=str(int(float(value)))
				cv = CVTerm();
				cv.setQualifierType(BIOLOGICAL_QUALIFIER);
				cv.setBiologicalQualifierType(BQB_IS);
				cv.addResource(uri_compounds[key]+value);
				s1.addCVTerm(cv);




def create_parameters():
	for key in parameters.keys():
		value=parameters[key]
		parameter=model.createParameter()
		parameter.setId(key)
		parameter.setValue(value)
		parameter.setConstant(True)


def create_reactions():
	for rxnid, row in rxns.iterrows():
		print(rxnid)
		# need to place in notes why reaction included
		if rxnid not in inclusion_rxns:
			continue
		rxnid=str(rxnid)
		rxn_name=str(rxns['Reaction'][rxnid])
		subsystem=rxns['Subsystem'][rxnid]
		if 'exchange reaction' in subsystem:
			subsystem='Exchange reaction'
		reversible=rxns['Reversible'][rxnid]
		ub=rxns['UB'][rxnid]
		lb=rxns['LB'][rxnid]
		objective=rxns['Objective'][rxnid]
		sbo=str(rxns['SBO'][rxnid])
		# literature
		notes=rxns['Notes'][rxnid]
		references=str(rxns['References'][rxnid])
		# reaction links
		ec_code=str(rxns['E.C. number'][rxnid])
		rhea=str(rxns['rhea'][rxnid])
		metacyc_reaction=str(rxns['metacyc'][rxnid])
		kegg_reaction=str(rxns['kegg'][rxnid])
		bigg_reaction=str(rxns['bigg'][rxnid])
		mnx_reaction=str(rxns['MetaNetX'][rxnid])
		# gene annotations
		complexes=rxns['Complex'][rxnid]
		fog_annotation=rxns['FOG'][rxnid]
		bypass_organism=rxns['organism'][rxnid]
		bypass_model=rxns['model'][rxnid]
		#
		print('basic information')
		r1 = model.createReaction()
		status=r1.setId('R_'+rxnid)
		check_status(status,'reaction_id','R_'+rxnid)
		status=r1.setMetaId('R_'+rxnid)
		check_status(status,'reaction_metid','R_'+rxnid)
		status=r1.setName(rxn_name)
		check_status(status,'reaction_name',rxn_name)
		status=r1.setSBOTerm(sbo)
		if reversible==1:
			r1.setReversible(True)
		elif reversible==0:
			r1.setReversible(False)
		r1.setFast(False)
		notes_list=['<p>SUBSYSTEM: '+subsystem+'</p>']
		if len(notes)>0:
			notes_list.append('<p>NOTES: '+str(notes.replace('\n','$').replace('[','_').replace(']','_').replace('?','').replace('<=>','rxn_rev'))+'</p>')
		if len(bypass_model)>0:
			notes_list.append('<p>MODEL_OVERRIDE: '+bypass_model+'</p>')
		if len(bypass_organism)>0:
			notes_list.append('<p>ORGANISM_OVERRIDE: '+bypass_organism+'</p>')
		r1.appendNotes('<notes><body xmlns="http://www.w3.org/1999/xhtml">'+str(''.join(notes_list))+'</body></notes>')
		#need to add in fugure
		#if oid_rxns[oid][rxnid] not in ['FOGS']:
		#	r1.appendNotes('<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>'+str(oid_rxns[oid][rxnid])+'</p></body></notes>')
		rplugin = r1.getPlugin("fbc")
		if rplugin != None:
			rplugin.setLowerFluxBound( [key for key in parameters if parameters[key]==lb][0] )
			rplugin.setUpperFluxBound( [key for key in parameters if parameters[key]==ub][0] )
		print('database links')
		reaction_identifiers=[
			('rhea',rhea),
			('kegg.reaction',kegg_reaction),
			('metacyc.reaction',metacyc_reaction),
			('ec-code',ec_code),
			('mnx.reaction',mnx_reaction),
			('bigg.reaction',bigg_reaction)]
		# identifiers
		for key,values in sorted(reaction_identifiers):
			if len(values)>0:
				values=[str(value.split(':')[0]) for value in values.split('|')]
				for value in values:
					cv = CVTerm();
					cv.setQualifierType(BIOLOGICAL_QUALIFIER);
					cv.setBiologicalQualifierType(BQB_IS);
					cv.addResource(uri_reactions[key]+value);
					r1.addCVTerm(cv);
		# literature references
		print('literature links')
		for reference in filter(None,references.split(';')):
			key=reference.split(':')[0]
			values=reference.split(':')[1].split('|')
			for value in values:
				cv = CVTerm();
				cv.setQualifierType(BIOLOGICAL_QUALIFIER);
				cv.setBiologicalQualifierType(BQB_IS_DESCRIBED_BY);
				cv.addResource(uri_literature[key]+value);
				r1.addCVTerm(cv);
		# parse reaction
		print('parse reactions')
		rxn_formula=rxns['Formula'][rxnid]
		parsed_rxn=parse_formula(rxn_formula)
		# reactants
		for stoich,met in parsed_rxn['reactants']:
			reactant = r1.createReactant()
			reactant.setSpecies('M_'+met[:-3]+'_'+met[-2:-1])
			reactant.setStoichiometry(float(stoich))
			reactant.setConstant(True)
		# products
		for stoich,met in parsed_rxn['products']:
			print('product')
			product = r1.createProduct()
			product.setSpecies('M_'+met[:-3]+'_'+met[-2:-1])
			product.setStoichiometry(float(stoich))
			product.setConstant(True)
		# genes
		# this section can use some major clean-up
		print('parse genes')
		if len(fog_annotation)>0:
			#gene_associations=process_fog_annotation(oid,fog_annotation)
			complex_variants=[';'.join(complex_variant) for oid in inclusion_oids for complex_variant in process_fog_annotation(oid,fog_annotation).values()]
			complex_variants=set(complex_variants)
			# need to account for complex annotations that are forced via organisms inclusion
			bypass_list=['BIOMASS','DEMAND','DIFFUSION','EQUILIBRIUM','ESSENTIAL','EXCHANGE','SPONTANEOUS','ORPHAN']
			if len(complex_variants)==0 and ( any(oid in inclusion_oids for oid in bypass_organism.split('|')) or bypass_model in bypass_list):
				complex_variants=[ ';'.join([subunit for subcomplex in variant.split('&&') for subunit in subcomplex.split('&')]) for variant in fog_annotation.replace('(','').replace(')','').replace(' ','').split('||') ]
				complex_variants=[';'.join([fog for fog in variant.split(';') if any( len(aybrah[oid][fog[:8]])>0  for oid in inclusion_oids)     ])  for variant in complex_variants]
				complex_variants=filter(None,complex_variants)
			elif len(complex_variants)>0:
				# sort complex variants by size
				_,complex_variants=zip(*sorted([(cv.count(';'),cv) for cv in complex_variants]))
			# test is not complex annotaiton
			isozyme_fogs=[complex_variant for complex_variant in complex_variants if ';' not in complex_variant]
			# isozyme annotation, ie not a complex annotation
			if len(isozyme_fogs)==len(list(complex_variants)) and len(isozyme_fogs)>1:
				rplugin = r1.getPlugin("fbc")
				gpr=rplugin.createGeneProductAssociation()
				g_current_level1=gpr.createOr()
				for fog in isozyme_fogs:
					# if not strain just use FOG
					if level!='Strain':
						gene_tag=g_current_level1.createGeneProductRef()
						gene_tag.setGeneProduct(str(fog[:8]))
					# if strain populate ORF
					else:
						genes=[gene.split('|')[1] for gene in aybrah[oid][fog[:8]].split(';')]
						for gene in genes:
							gene_tag=g_current_level1.createGeneProductRef()
							gene_tag.setGeneProduct('G_'+gene)
			elif len(list(complex_variants))>0:
				_,complex_variants=zip(*sorted([(complex_variant.count(';'),complex_variant)for complex_variant in complex_variants]))
				#gpr_instance_list=sorted(gene_associations.keys())
				rplugin = r1.getPlugin("fbc")
				gpr=rplugin.createGeneProductAssociation()
				# add complex in notes
				cpxs=complexes.replace('&','||').replace('(','').replace(')','').replace(' | ','||').replace(' ','').split('||')
				for c in set(cpxs):
					if len(c)==0:
						continue
					if c[:3]=='CPX':
						note = '<body xmlns="http://www.w3.org/1999/xhtml"><p>'+'http://identifiers.org/complexportal/'+str(c)+'</p></body>'
					else:
						note = '<body xmlns="http://www.w3.org/1999/xhtml"><p>'+'AYbRAH complex: '+str(c)+'</p></body>'
					gpr.appendNotes(str(note))
					# cannot get CV term added to geneProductAssociation
					#for c in complexes.replace('&','||').replace('(','').replace(')','').replace(' ','').split('||'):
					#	cv = CVTerm();
					#	cv.setQualifierType(BIOLOGICAL_QUALIFIER);
					#	cv.setBiologicalQualifierType(BQB_IS_DESCRIBED_BY);
					#	cv.addResource('http://identifiers.org/complexportal/'+str(c));
					#	gpr.addCVTerm(cv);
				if len(complex_variants)>1:
					# multiple instances of complexes/enzymes
					g_current_level1=gpr.createOr()
				else:
					g_current_level1=gpr
				for complex_variant in complex_variants:
					fogs=complex_variant.split(';')
					if len(fogs)>1:
						# tag needed for complex
						g_current_level2=g_current_level1.createAnd()
					else:
						# just one gene
						g_current_level2=g_current_level1
					for fog in fogs:
						fog
						# if not strain just use FOG
						if level!='Strain':
							gene_tag=g_current_level2.createGeneProductRef()
							gene_tag.setGeneProduct(str(fog[:8]))
							if '?' in fog:
								gene_tag.setAnnotation('non-essential')
						# if strain populate ORF
						else:
							genes=[gene.split('|')[1] for gene in aybrah[oid][fog[:8]].split(';')]
							if len(genes)==1:
								gene_tag=g_current_level2.createGeneProductRef()
								gene_tag.setGeneProduct('G_'+genes[0])
							else:
								g_current_level3=g_current_level2.createOr()
								for gene in genes:
									gene_tag=g_current_level3.createGeneProductRef()
									gene_tag.setGeneProduct('G_'+gene)






def create_objective():
	mplugin = model.getPlugin("fbc");
	#
	objective = mplugin.createObjective();
	objective.setId("objective");
	objective.setType("maximize");
	#
	mplugin.setActiveObjectiveId("objective");
	fluxObjective = objective.createFluxObjective();
	fluxObjective.setReaction("R_BIOMASS");
	fluxObjective.setCoefficient(1);


def check_status(status,tag,entry):
	if status!=0:
		fname=taxonomy['name'][index]+'.txt'
		open(fname,'a').write('\t'.join([str(status),tag,entry])+'\n')



def create_geneproducts(uniprot_proteome):
	print('parse gene reactions')
	# this is adding all of them
	fogs=filter(None,'|'.join(rxns['FOG'].tolist()).replace('(','').replace(')','').replace(' ','').replace('&','|').split('|'))
	fogs=sorted(list(set([fog[:8] for fog in fogs])))
	mplugin = model.getPlugin("fbc");
	# load uniprot, could also use strain variable but not all strains are populated
	#if len(uniprot_proteome)>0:
	#	tree=ET.parse(urllib2.urlopen('https://www.uniprot.org/uniprot/?query=proteome:'+uniprot_proteome+'&format=xml'))
	#	root = tree.getroot()
	#else:
	#	root={}
	#
	for fog in fogs:
		#print fog
		fog=str(fog[:8])
		if fog in fog_descriptions.index:
			protein_id=str(fog_descriptions['Protein ID'][fog])
			protein_name=str(fog_descriptions['Protein description'][fog])
		else:
			protein_id,protein_name='',''
		#print '\t'+protein_name
		if level=='Strain':
			seqids=list(filter(None,[seqid for seqid in aybrah[inclusion_oids[0]][fog].split(';')]))
			for seqid in seqids:
				print(fog,seqid)
				locus_tag=seqid.split('|')[1]
				gene = mplugin.createGeneProduct();
				status=gene.setId('G_'+locus_tag);
				check_status(status,'gene_id',locus_tag)
				gene.setSBOTerm('SBO:0000243')
				#status=gene.setMetaId('G_'+locus_tag);
				#check_status(status,'gene_metaid',locus_tag)
				if len(protein_name)>0:
					gene.setName(protein_name)
					check_status(status,'gene_name',protein_name)
				# need to check if protein ID exists (such as ALD6.3).
				# if there are paralogs need to append character to make unique
				if len(protein_id)>0 and len(seqids)>1:
					letter=chr(65+seqids.index(seqid))
					gene.setLabel(protein_id+'.'+letter)
					check_status(status,'gene_label',protein_id+'.'+letter)
				elif len(protein_id)>0 and len(seqids)==1:
					status=gene.setLabel(protein_id)
					check_status(status,'gene_label',protein_id)
				else:
					status=gene.setLabel('G_'+locus_tag);
					check_status(status,'gene_label','G_'+locus_tag)
				#
				#accessions_kegg,accessions_refseq,accessions_gids=get_attributes(locus_tag)
				#for value in accessions_kegg:
				#	cv = CVTerm();
				#	cv.setQualifierType(BIOLOGICAL_QUALIFIER);
				#	cv.setBiologicalQualifierType(BQB_IS_ENCODED_BY);
				#	url=uri_genes['kegg.genes']+value
				#	cv.addResource(url.replace('https','http'));
				#	gene.addCVTerm(cv);
		else:
			seqids=filter(None,[seqid for oid in inclusion_oids for seqid in aybrah[oid][fog].split(';')])
			if len(seqids)>0:
				gene = mplugin.createGeneProduct();
				gene.setId(fog);
				#gene.setMetaId(fog);
				if len(protein_name)>0:
					gene.setName(protein_name)
				if len(protein_id)>0:
					gene.setLabel(protein_id)
				else:
					gene.setLabel(fog)







def get_attributes(locus_tag):
	for entry in root.getchildren():
		for accession in entry.findall('{http://uniprot.org/uniprot}accession'):
			if accession.text == locus_tag:	
				accessions_kegg=[db.attrib['id'] for db in entry.findall('{http://uniprot.org/uniprot}dbReference') if db.attrib['type']=='KEGG']
				accessions_refseq=[db.attrib['id'] for db in entry.findall('{http://uniprot.org/uniprot}dbReference') if db.attrib['type']=='RefSeq']
				accessions_gids=[db.attrib['id'] for db in entry.findall('{http://uniprot.org/uniprot}dbReference') if db.attrib['type']=='GeneID']
				return accessions_kegg,accessions_refseq,accessions_gids






# what about gfp notes and others?


# other
uri_other={
	'bigg_compartment':'http://identifiers.org/bigg.compartment/',
	'biocyc':'http://identifiers.org/biocyc/'}




uri_genes={
	'refseq':'https://identifiers.org/refseq',
	'uniprot':'https://identifiers.org/uniprot',
	'kegg.genes':'https://identifiers.org/kegg.genes/',
	'ncbigi':'https://identifiers.org/ncbigi/',
	'ncbigene':'https://identifiers.org/ncbigene/',
	'ncbiprotein':'https://identifiers.org/ncbiprotein/'
}



# reaction
uri_reactions={
	'rhea':'http://identifiers.org/rhea/', # done
	'metacyc.reaction':'http://identifiers.org/metacyc.reaction/', # done
	'kegg.reaction':'http://identifiers.org/kegg.reaction/', # done
	'bigg.reaction':'http://identifiers.org/bigg.reaction/', # done
	'ec-code':'http://identifiers.org/ec-code/', # done
	'mnx.reaction':'http://identifiers.org/metanetx.reaction/'} # done

# compound
uri_compounds={
	'chebi':'http://identifiers.org/chebi/CHEBI:',
	'metacyc.compound':'http://identifiers.org/metacyc.compound/',
	'bigg.metabolite':'http://identifiers.org/bigg.metabolite/',
	'kegg.compound':'http://identifiers.org/kegg.compound/',
	'pubchem.compound':'http://identifiers.org/pubchem.compound/',
	'inchi':'http://identifiers.org/inchi/',
	'inchikey':'http://identifiers.org/inchikey/',
	'mnx.species':'http://identifiers.org/metanetx.chemical/',
	'bigg.species':'http://identifiers.org/bigg.metabolite/'}


uri_compartments={
	'mnx_compartment':'http://identifiers.org/metanetx.compartment/',
	'bigg_compartment':'http://identifiers.org/bigg.compartment/'
}

# literature
uri_literature={
	'doi':'http://identifiers.org/doi/',
	'uspto':'http://identifiers.org/uspto/',
	'google_patent':'http://identifiers.org/google.patent/',
	'isbn':'http://identifiers.org/isbn/',
	'pubmed':'http://identifiers.org/pubmed/'}



"""

complex information
'https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000657'


ChEBI
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000002


MetaCyc Compound
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000609

MetaCyc Reaction
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000610

BIGG Compartment
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000555

BIGG reaction
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000557

BIGG metabolite
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000556

KEGG Compound
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000013

Google patent
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000537

ISBN
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000064

USPTO
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000538

DOI
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000019

RHEA
https://www.ebi.ac.uk/miriam/main/datatypes/MIR:00000082

PubChem
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000034

BioCyc
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000194

PubMed
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000015

InChI
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000383

KEGG Reaction
https://www.ebi.ac.uk/miriam/main/collections/MIR:00000014

"""

root_aybrah='./aybrah'

# local_path
#root_aybrah='https://github.com/kcorreia/aybrah/raw/master/'


path_aybrah_xlsx=root_aybrah+'/aybrah.xlsx'


taxonomy=pd.read_excel(path_aybrah_xlsx,sheet_name='taxon_nodes').fillna('')
#taxonomy.drop(46)
drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)

_,taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))


version_aybrah=open(root_aybrah+'/version.txt','r').read()
aybrah=pd.read_csv(root_aybrah+'/aybrah.tsv',sep='\t').fillna('').set_index('FOG')


organisms=pd.read_excel(path_aybrah_xlsx,sheet_name='curated_taxonomy_fungi').fillna('')






"""
Should have quality checks for all lb/ub, no duplicates in rxnid
no blank mets
"""





root_fyrment='./fyrment'
#url_fyrment_xlsx='https://github.com/LMSE/fyrment/raw/master/fyrment.xlsx'

version_fyrment=open(root_fyrment+'/version.txt','r').read()

path_fyrment_xlsx=root_fyrment+'/fyrment.xlsx'
#encoding='utf-8'

# used to add descriptions for gene products in SBML
# fbc:label and fbc:name
# <fbc:geneProduct sboTerm="SBO:0000243" fbc:id="G_Q6CWR9" fbc:name="RPL1 Large ribosome protein subunit" fbc:label="RPL1"/>
fog_descriptions=pd.read_csv(root_aybrah+'/descriptions_fog.txt',sep='\t').set_index('FOG').fillna('')


# the FYRMENT in Excel format has 3 sheets:
# rxns that include FOG annotations
rxns=pd.read_excel(path_fyrment_xlsx,sheet_name='reactions',skiprows=[0]).fillna('').set_index('Rxn name')
# metabolites that include species information like formula, charge, meta data
mets=pd.read_excel(path_fyrment_xlsx,sheet_name='metabolites').fillna('').set_index('Metabolite name')
# compartment information
df_c=pd.read_excel(path_fyrment_xlsx,sheet_name='compartments').set_index('compartment_id').fillna('')

# organism IDs in AYbRAH
oids=aybrah.columns[2:].tolist()


sbml_level=3
sbml_version=1

parameters={
	'cobra_default_inf':1000,
	'cobra_default_neginf':-1000,
	'cobra_default_substrate':-10,
	'cobra_default_zero':0
}



# creates a dictionary of which rxnids to include for each oid in the network
oid_rxns=get_rxn_inclusion()




# create models for all taxonomies (from kingdom to strain)
for index,row in taxonomy.iterrows():
# create models for all strains
#for index,row in taxonomy[taxonomy['level']=='Strain'].iterrows():
# index 113 for K. lactis metabolic model
#for index in [113]:
	#break
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	name=taxonomy['name'][index]
	readable=taxonomy['readable'][index].replace('.','_')
	#
	path_xml_dir = './genre/xml/'+str(level_order).zfill(2)+'_'+level.lower()
	path_xml_file=path_xml_dir+'/i'+readable+'.xml'
	if os.path.isfile(path_xml_file):
		print('skip '+readable)
		continue
	else:
		print('build '+readable)
	#name
	inclusion_oids=taxonomy['oids'][index].split('|')
	# only include oid annotated in aybrah
	inclusion_oids=[oid for oid in inclusion_oids if oid in aybrah.columns[2:]]
	if name=='Eukaryota' or name=='Fungi':
		inclusion_rxns=rxns.index.tolist()
	elif len(inclusion_oids)==0:
		continue
	else:
		inclusion_rxns=set([rxnid for oid in inclusion_oids for rxnid in oid_rxns[oid].keys() if oid in oid_rxns.keys()])
	# used to debug stoich consistency
	#inclusion_rxns=[r for r in inclusion_rxns if r not in skip]
	#
	model_id=name.replace(' ','_').replace('-','_')+'_GENRE_FYRMENT_'+version_fyrment.replace('.','_')
	if level=='Strain':
		#continue
		model_metaid=name.replace(' ','_')+'_GENRE_FYRMENT_'+version_fyrment
		model_name=name+' genome-scale network reconstruction'
		oid=taxonomy['oids'][index]
		uniprot_proteome=taxonomy['uniprot:proteome'][index]
	else:
		#continue
		model_metaid=name.replace(' ','_')+'_GENRE_FYRMENT_'+version_fyrment
		model_name=name+' genome-scale network reconstruction'
		oid=None
		uniprot_proteome=None
	#if os.path.exists(model_id+'.xml'):
	#	continue
	###################################
	# header information
	sbmlns = SBMLNamespaces(3,1)
	sbmlns.addPackageNamespace('fbc',2)
	sbmlns.addPackageNamespace('groups',1)
	document = SBMLDocument(sbmlns)
	document.setPackageRequired("fbc", False)
	document.setPackageRequired("groups", False)
	namespaces = document.getNamespaces()
	status = namespaces.add('http://www.sbml.org/sbml/level3/version1/fbc/version2', 'fbc')
	status = namespaces.add('http://www.sbml.org/sbml/level3/version1/groups/version1', 'groups')
	model = document.createModel()
	mplugin = model.getPlugin("fbc");
	mplugin.setStrict(True)
	# Set model name
	#
	model.setId(str(model_id))
	model.setMetaId(str(model_metaid))
	model.setName(str(model_name))
	# include AYbRAH, FYRMENT versions
	notes=['<p>'+db+' DB version: '+str(v)+'</p>' for db,v in zip(['AYbRAH','FYRMENT'],[version_aybrah,version_fyrment])]
	note = '<body xmlns="http://www.w3.org/1999/xhtml">'+''.join(notes)+'</body>'
	model.appendNotes(note)
	#
	###################################
	create_units()
	###################################
	# compartment information
	# based on inclusion_rxns
	metabolites_in_reactions=get_metabolites_in_reactions()
	compartments=set([met[-2:-1] for met in set(metabolites_in_reactions)])
	"""
	Test for current metabolites
	for a in set([(met[-2:-1],met) for met in set(output_reaction_metabolites)]):
		open('compartments2.txt','a').write('\t'.join(a)+'\n')
	"""
	create_compartments()
	###################################
	# metabolite species
	create_species(name)
	###################################
	# reactions
	create_reactions()
	#writeSBML(document,model_id+'.xml')
	###################################
	#
	create_parameters()
	#
	###################################
	#
	create_objective()
	#
	###################################
	#
	create_geneproducts(uniprot_proteome)
	#
	###################################
	if not os.path.exists(path_xml_dir):
		os.makedirs(path_xml_dir)
	#raise Exception ('review before saving')
	#path_xml_file='test_kla.xml'
	writeSBML(document,str(path_xml_file))
	###################################




#df_dna=pd.read_csv('./biomass/dna_stoich.txt',sep='\t').set_index('oid')
df_dna=pd.read_excel('./biomass/dna_stoich.xlsx').set_index('oid')


# update DNA composition, and other organism specific constraints
# can/should update this function such that it parses external files rather than hard coding changes


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
	elif level=='Strain':
		print(name)
	else:
		continue
	#
	reader = SBMLReader()
	#
	document = reader.readSBML(str(path_xml_file))
	model = document.getModel()
	rxnids_all=[reaction.id for reaction in model.getListOfReactions()]
	#
	for reaction in model.getListOfReactions():
		if reaction.id=='R_BIOMASS_DNA':
			# break
			# update the biomass based on GC content
			for reactant in reaction.reactants:
				metid=reactant.species[2:-2]
				stoich=df_dna[metid][oid]
				print(metid,stoich)
				#reactant.stoichiometry=str(abs(df_dna[metid][oid]))
				reactant.setStoichiometry(abs(df_dna[metid][oid]))
			for product in reaction.products:
				metid=product.species[2:-2]
				stoich=df_dna[metid][oid]
				print(metid,stoich)
				#reactant.stoichiometry=str(df_dna[metid][oid])
				product.setStoichiometry(df_dna[metid][oid])
		# if inositol auxtrophy exists, turn on inositol exchanage
		if reaction.id=='R_EX_inost_e' and all(x not in rxnids_all for x in ['R_MI1PS','R_MI1PP']):
			print('inositol auxotrophy')
			rplugin = reaction.getPlugin("fbc")
			rplugin.setLowerFluxBound('cobra_default_neginf')
		# if inositol degradation pathway exists and auxotrophic, shut off degradation pathway
		if reaction.id=='R_INOSTO' and all(x not in rxnids_all for x in ['R_MI1PS','R_MI1PP']):
				print('turn off inositol degradation pathway')
				rplugin = reaction.getPlugin("fbc")
				rplugin.setLowerFluxBound('cobra_default_zero')
				rplugin.setUpperFluxBound('cobra_default_zero')	
			# "mitochondrial origin of aspartate biosynthesis in Yarrowia lipolytica"
			# https://academic.oup.com/femsyr/article/5/6-7/545/528244
		if reaction.id=='R_ASPTA' and oid == 'yli':
			print('turn off cytosolic aspartate biosynthesis')
			rplugin = reaction.getPlugin("fbc")
			rplugin.setLowerFluxBound('cobra_default_zero')
			rplugin.setUpperFluxBound('cobra_default_zero')
		# alternative methionine metabolism in Ogataea
		# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100725
		if reaction.id=='R_CYSTL' and oid=='opm':
			print('turn on met__L biosynthesis from cystL')
			rplugin = reaction.getPlugin("fbc")
			rplugin.setLowerFluxBound('cobra_default_zero')
			rplugin.setUpperFluxBound('cobra_default_inf')
		# "cytosolic origin of alanine biosynthesis in S. kluyveri"
		# https://academic.oup.com/femsyr/article/5/6-7/545/528244
		# need inositol pathways shut off for inositol auxotrophic strains
	writeSBML(document,str(path_xml_file))















