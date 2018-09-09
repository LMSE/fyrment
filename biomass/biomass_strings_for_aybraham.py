import pandas as pd


"""

AMINO ACIDS (aa)

basic peptide / protein equation:
sum_of(v_aa * AA, for all aa) + net_charge h+   -> 1 protein (of some MW) + sum_of(v_aa,for all aa)-1 h2o

we can adjust the v_aa, while keeping the composition constant, such that:
-protein net charge is 0
-protein MW is 1000 g/mmol

Note A: coefficient does not apply to moles_2 of water
Example 3 ala -> 2 h2o + peptide (3-ala)
		6 ala -> 5 h2o + peptide (6-ala)

Note B:
net_charge = sum_of(v_aa*Q_aa)
MW protein = sum_of(v_aa*MW_aa)-net_charge -[sum_of(v_aa)-1]*MW_h2o

Note C:
We want to get the scaling factor (beta) such that MW protein = 1000 g/mmol:

analytical solution (see corresponding Excel spreadsheet):
alpha = sum_of(v_aa*MW_aa) - net_charge - sum_of(v_aa)*M_h2o
beta = 1000/MW_1 (  1 + MW_h2o / alpha ) - MW_h2o / alpha

Note D:
2 GTP are required for amino acid
Therefore net stoich for tRNA equation
GTP: 2*moles_2 (reactant)
H2O: 2*moles_2-(moles_2-1) = moles_2+1 (reactant)
GDP: 2*moles_2 (product)
P:   2*moles_2 (product)
H:   2*moles_2+net_charge_2 (product)

"""

def protein_formulas():
	"""
	Returns reaction formula with simple amino acids as reactants and with tRNA + 2 GTP per amino acid.
	Input data is mol amino acid / mol protein.
	"""
	df_calc = pd.DataFrame({'MW':aa_table['mass'],'charge':aa_table['charge']})
	net_charge_1=sum([aa_table['iMM904 mmol AA/gDCW'][aa]*aa_table['charge'][aa] for aa,row in aa_table.iterrows()])
	MW_1=sum([aa_table['iMM904 mmol AA/gDCW'][aa]*aa_table['mass'][aa] for aa,row in aa_table.iterrows()])\
			-net_charge_1*1.01-( sum(aa_table['iMM904 mmol AA/gDCW'].tolist()) -1)*18.02
	# See Note C
	alpha=sum([aa_table['iMM904 mmol AA/gDCW'][aa]*aa_table['mass'][aa] for aa,row in aa_table.iterrows()])\
			-net_charge_1*1.01-sum([aa_table['iMM904 mmol AA/gDCW'][aa] for aa,row in aa_table.iterrows()])*18.02
	beta=1000.0/MW_1*(1+18.02/alpha)-18.02/alpha
	#
	df_calc['mmol/gDCW']=[-beta*v_aa for v_aa in aa_table['iMM904 mmol AA/gDCW'].tolist()]
	moles_2=-sum(df_calc['mmol/gDCW'].tolist())
	net_charge_2=sum([df_calc['mmol/gDCW'][aa]*df_calc['charge'][aa] for aa,row in df_calc.iterrows()])
	# species for simple protein equation with GTP cost
	df_simple=df_calc.copy()
	df_simple.loc['h']=[1.01,1,net_charge_1*beta]
	df_simple.loc['h2o']=[18.02,0,-sum(df_calc['mmol/gDCW'].tolist())-1]
	MW_2=round(-sum([df_simple['MW'][species]*df_simple['mmol/gDCW'][species] for species,row in df_simple.iterrows()]),4)
	df_simple.loc['protein_1g']=[MW_2,0,1]
	df_trna=df_calc.copy()
	df_trna['charged']=aa_table['charged'].tolist()
	#df_trna=df_trna.reset_index()
	df_trna=df_trna.set_index('charged')
	# species for tRNA-based protein equation with GTP cost
	df_trna.loc['gtp']=[519.14886,-4,-(2*moles_2)]
	df_trna.loc['h2o']=[18.02,0,-(2*moles_2+1)]
	for aa_i,row in df_calc.iterrows():
		#aa_table['uncharged'][aa_i]
		aa_uncharged=aa_table['uncharged'][aa_i]
		aa_charged=aa_table['charged'][aa_i]
		df_trna.loc[aa_uncharged]=['','',-df_trna['mmol/gDCW'][aa_charged]]
	df_trna.loc['gdp']=[440.17670,-3,(2*moles_2)]
	df_trna.loc['pi']=[95.9793,-2,(2*moles_2)]
	df_trna.loc['h']=[1.01,1,3*moles_2-net_charge_2]
	df_trna.loc['protein_1g']=[MW_2,0,1]
	df_simple.to_csv('output_table_protein_simple.txt',sep='\t')
	df_trna.to_csv('output_table_protein_trna.txt',sep='\t')
	return df_simple,df_trna


# amino acid composition from several sources
aa_table=pd.read_csv('data_aa_composition.txt',sep='\t',comment='#').set_index('amino acid').fillna('0')


df_simple,df_trna=protein_formulas()

for df in [df_simple,df_trna]:
	for compartment in list('cm'):
		reactants=' + '.join([str(-round(df['mmol/gDCW'][species],4))+' '+species+'['+compartment+']' for species,row in df[df['mmol/gDCW']<0].iterrows()])
		products=' + '.join([str(round(df['mmol/gDCW'][species],4))+' '+species+'['+compartment+']' for species,row in df[df['mmol/gDCW']>0].iterrows()])
		formula=reactants+' -> '+products
		open('formula_protein.txt','a').write('\t'.join(['Protein synthesis ['+compartment+']',formula])+'\n')



###########################################################################################
#
#	RNA,DNA
#	i=ATP,GTP,CTP,UTP or dATP,dGTP,dCTP,dUTP
#	Mass balance of RNA / DNA polymerase reaction
#	sum(vi*MW_i) + sum(vi)*MW_H + MW_H2O = MW_RNA + sum(vi)*MW_ppi
#	adj = value used to mulplipy vi such that MW-RNA = 1000 mg/mmol
#	1 mol H2O required to remove ppi from 5' end of RNA/DNA	
#
##174.95126 ppi
###########################################################################################
#
#
#
#





xtp=pd.read_csv('data_xtp_composition.txt',sep='\t').fillna('').set_index('metid_reactant')


data_sources=xtp.columns[xtp.columns.tolist().index('PIC RNA'):].tolist()


def XTP_formulas_from_table(macromolecule,source):
	df_temp=xtp[xtp['macro']==macromolecule].copy()
	# get only rows with values
	#df_temp=xtp[xtp[data_source].apply(lambda x: isinstance(x, float))]
	#macromolecule=df_temp['macro'][0]
	n=sum(df_temp[source])
	MW_nucleotides_delta1=sum([v*MW for v,MW in zip(df_temp[source],df_temp['MW_reactant'])])+\
		sum(df_temp[source])*1.01-sum(df_temp[source])*174.95126
	# see RNA comment section to see how adj is derived
	adj=(1000-18.02)/(MW_nucleotides_delta1)
	df_calc = pd.DataFrame({'MW':df_temp['MW_reactant'],'charge':df_temp['charge_reactant']})
	df_calc['mmol/gDCW']= [round(-adj*v,3) for v in df_temp[source]]
	df_calc.loc['h']=[1.01,1,round(-adj*n,3)]
	df_calc.loc['h2o']=[18.02,0,-1]
	df_calc.loc['ppi']=[174.95126,-3,round(adj*n,3)]
	df_calc.loc[macromolecule.lower()+'_1g']=[1000,0,1]
	check=sum([v*MW for v,MW in zip(df_calc['mmol/gDCW'],df_calc['MW'])])
	df_calc.to_csv('output_table_'+macromolecule+'_'+source+'.txt',sep='\t')
	# error on 1 gram
	print('Error '+str(round(check/1000*100,4))+'%')
	return df_calc




df_rna=XTP_formulas_from_table('RNA','iMM904 mmol XTP/gDCW')
df_dna=XTP_formulas_from_table('DNA','iMM904 mmol XTP/gDCW')

df_rna=XTP_formulas_from_table('RNA','Sanchez2017')
df_dna=XTP_formulas_from_table('DNA','Sanchez2017')




for df in [df_rna,df_dna]:
	for compartment in list('cm'):
		reactants=' + '.join([str(-round(df['mmol/gDCW'][species],4))+' '+species+'['+compartment+']' for species,row in df[df['mmol/gDCW']<0].iterrows()])
		products=' + '.join([str(round(df['mmol/gDCW'][species],4))+' '+species+'['+compartment+']' for species,row in df[df['mmol/gDCW']>0].iterrows()])
		formula=reactants+' -> '+products
		open('formula_xtp.txt','a').write('\t'.join(['XNA synthesis ['+compartment+']',formula])+'\n')



# formula used to calculate dXTP for specific genomes

def dXTP_coefficients_from_GC_content(GC_content):
	macromolecule='DNA'
	df_temp=xtp[xtp['macro']==macromolecule].copy()
	df_temp['genomic']=''
	df_temp['genomic']['datp']=(1-GC_content)/2
	df_temp['genomic']['dttp']=(1-GC_content)/2
	df_temp['genomic']['dgtp']=GC_content/2
	df_temp['genomic']['dctp']=GC_content/2
	# get only rows with values
	#df_temp=xtp[xtp[data_source].apply(lambda x: isinstance(x, float))]
	#macromolecule=df_temp['macro'][0]
	n=sum(df_temp['genomic'])
	MW_nucleotides_delta1=sum([v*MW for v,MW in zip(df_temp['genomic'],df_temp['MW_reactant'])])+\
		sum(df_temp['genomic'])*1.01-sum(df_temp['genomic'])*174.95126
	# see RNA comment section to see how adj is derived
	adj=(1000-18.02)/(MW_nucleotides_delta1)
	df_calc = pd.DataFrame({'MW':df_temp['MW_reactant'],'charge':df_temp['charge_reactant']})
	df_calc['mmol/gDCW']= [round(-adj*v,3) for v in df_temp['genomic']]
	df_calc.loc['h']=[1.01,1,round(-adj*n,3)]
	df_calc.loc['h2o']=[18.02,0,-1]
	df_calc.loc['ppi']=[174.95126,-3,round(adj*n,3)]
	df_calc.loc[macromolecule.lower()+'_1g']=[1000,0,1]
	check=sum([v*MW for v,MW in zip(df_calc['mmol/gDCW'],df_calc['MW'])])
	df_calc.to_csv('output_table_genomic_'+macromolecule+'.txt',sep='\t')
	# error on 1 gram
	print('Error '+str(round(check/1000*100,4))+'%')
	return df_calc


url_aybrah_xlsx='https://github.com/kcorreia/aybrah/raw/master/aybrah.xlsx'
df_genomes=pd.read_excel(url_aybrah_xlsx,sheet_name='organisms').fillna('').set_index('order')



df_dna_reaction=pd.DataFrame(columns='datp|dttp|dctp|dgtp|h|h2o|ppi|dna_1g'.split('|'))

# how was GC content extracted? assembly accessions via UniProt but unsure where code is
for order,row in df_genomes.iterrows():
	oid=df_genomes['oid'][order]
	GC_content=df_genomes['GC content'][order]
	df_calc=dXTP_coefficients_from_GC_content(GC_content)
	df_dna_reaction.loc[oid]=df_calc['mmol/gDCW'].tolist()

df_dna_reaction.to_csv('aybraham_dna_stoich.txt',sep='\t')


###########################################################################################
#
#	ACOA and PLIPIDS
#	plipids are neutral molecueles, MW of R is calculated from specified acoa
#	stoich in plipid is used to define following reaction:
#	stoich * plipid mmol/gDCW/h -> 1g lipid/gDCW/h
#	PC_MONOUNSAT_1G	1g pc, monounsaturated	1.36 pc[c] -> pc_monounsat_1g[c]
#
#	the product of the above reaction serves as the reactant to a composite reaction of other phospholipids
#	1 gDW phospholipid composite, monounsaturated	0.3946 pc_monounsat_1g[c] + 0.3372 pe_monounsat_1g[c] + 0.0 pg_monounsat_1g[c] + 0.0145 pa_SC_monounsat_1g[c] + 0.0425 ps_monounsat_1g[c] + 0.1853 ptd1ino_monounsat_1g[c] + 0.026 clpn_monounsat_1g[c] -> phospholipid_1g[c]
#
###########################################################################################
#
#
#
#



def acoa_plipids_formulas(acoa_type):
	print(acoa_type)
	df_calc_plipids=plipids.copy()
	df_temp=acoa[acoa[acoa_type]>0]
	coeff={}
	coeff['mmol']=[]
	coeff['mmol']=[w/MW*1000 for w,MW in zip(df_temp[acoa_type],df_temp['mass'])]
	# mol% used to represent 1 acoa
	coeff['mol%']=[round(c/sum(coeff['mmol']),4) for c in coeff['mmol']]
	MW_acoa=round(sum(    [v*M for v,M in zip(        coeff['mol%'],df_temp['acoa mass']  )  ]),4)
	# 763 is the mass of coa moiety, R is the side group added to the various phospholipid compounds
	MW_Ravg=round(MW_acoa-763.502-12.01-16+1.01,2)
	MW_plipids=[round(mass_ex_R+R*MW_Ravg,4) for mass_ex_R,R in zip(plipids['mass ex. R'].tolist(),plipids['R'])]
	df_calc_plipids['mass']=MW_plipids
	df_calc_plipids['stoich']=[1000/MW for MW in MW_plipids]
	return df_calc_plipids



acoa=pd.read_csv('data_acoa_composition.txt',sep='\t').set_index('acoa').fillna(0)

plipids=pd.read_csv('data_plipid_composition.txt',sep='\t').fillna(0).set_index('phospholipid')


acoa_types=['sat','monounsat','polyunsat']

output={acoa_type:acoa_plipids_formulas(acoa_type) for acoa_type in acoa_types}
acoa_plipids_formulas('monounsat')





