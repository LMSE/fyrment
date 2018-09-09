
from pandas import ExcelWriter
from openpyxl import load_workbook
import urllib
import pandas as pd

def process_fog_annotation(oid,fog_annotation):
	"""
	Could and should clean this up
	Purpose is to take in a fog_annotation in AYbRAHAM and oid
	Check if all essential subunits are present
	"""
	def process_complex_variants(oid,fog_annotation):
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
				print oid,'include complex '+str(i)
			else:
				# return empty
				print 'oid','not all subunits for complex '+str(i)
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





version_aybrah=urllib.urlopen('https://github.com/kcorreia/aybrah/raw/master/version.txt').read()
aybrah=pd.read_csv('https://github.com/kcorreia/aybrah/raw/master/aybrah.tsv',sep='\t').fillna('').set_index('FOG')


url_aybrah_xlsx='https://github.com/kcorreia/aybrah/raw/master/aybrah.xlsx'
organisms=pd.read_excel(url_aybrah_xlsx,sheet_name='curated_taxonomy_fungi').fillna('')
organisms=organisms.drop([6,21])


taxonomy=pd.read_csv('https://github.com/kcorreia/aybraham/raw/master/taxonomy_level_oids.txt',sep='\t')
drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)



# need to make sure ? fogs are captured

url_aybraham_xlsx='https://github.com/kcorreia/aybraham/raw/master/aybraham.xlsx'
rxns=pd.read_excel(url_aybraham_xlsx,sheet_name='reactions',encoding='ascii',skiprows=[0]).fillna('').set_index('Rxn name')


writer = ExcelWriter('./complexes/aybraham_complex_curation.xlsx')

fog_annotations=list(set(rxns[rxns.FOG.str.contains('&')].FOG.tolist()))
for number,fog_annotation in enumerate(fog_annotations):
	rxnids=';'.join(rxns[rxns.FOG==fog_annotation].index.tolist())
	formula=';'.join(rxns[rxns.FOG==fog_annotation].Formula.tolist())
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
	fogs=[fog[:8] for i in complex_check.keys() for fog in complex_check[i]]
	fogs_nonessential=[fog[:8] for fog in fogs if '?' in fog]
	columns=['Rxn name','Formula','Include']+[(i,str(fog)) for i in complex_check.keys() for fog in complex_check[i]]
	df_temp=pd.DataFrame(columns=columns)
	for order,row in organisms.iterrows():
		oid=organisms['oid'][order]
		complexes_parsed=process_fog_annotation(oid,fog_annotation)
		include=len(complexes_parsed)
		row=[rxnids,formula,include]+[0]*(len(columns)-3)
		df_temp.loc[oid]=row
		for num,fog in df_temp.columns[3:].tolist():
			fog
			seqids=aybrah[oid][fog[:8]]
			if len(seqids)>0:
				df_temp.loc[oid,(num,fog)]=1
	df_temp.to_excel(writer,str(number))


#df.loc[i, 'thumbnail'] = filename



writer.save()

"""
update colours to make it easier to identify missing subunits

"""

from openpyxl.styles import Color, PatternFill, Font, Border


yellowFill = PatternFill(start_color='FCFF00',
                   end_color='FCFF00',
                   fill_type='solid')

greyFill = PatternFill(start_color='a5a5a5',
                   end_color='a5a5a5',
                   fill_type='solid')




wb = load_workbook(filename = './complexes/aybraham_complex_curation.xlsx' )


for sheetname in wb.sheetnames:
	# highlight non-essential
	ws=wb[sheetname]
	for col in ws.iter_cols(min_row=1,max_row=1,min_col=2):
		for cell in col:
			if '?' in cell.value:
				cell.fill =yellowFill
	# grey out missing genes
	for col in ws.iter_cols(min_row=2):
		for cell in col:
			if str(cell.value)=='0':
				cell.fill = greyFill


wb.save('./complexes/aybraham_complex_curation.xlsx')


# https://stackoverflow.com/questions/30484220/python-fill-cells-with-colors-using-openpyxl
# ws['A1'].fill = redFill should work fine.

"""
	mismatching_rxn_columns=[cell.value \
			for header,col in zip(rxn_headers,ws.iter_cols(min_row=1, max_col=15, max_row=1)) \
			for cell in col if header != cell.value]
"""
