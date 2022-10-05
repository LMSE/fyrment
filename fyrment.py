


import pandas as pd


class FYRMENT:
	"""
	In development
	"""
	def __init__(self,path_fyrment_xlsx):
		print('\tLoad FYRMENT')
		#self.df=pd.read_excel('https://github.com/LMSE/fyrment/raw/master/fyrment.xlsx',sheet_name='reactions').fillna('').set_index('id')
		self.rxns=pd.read_excel(path_fyrment_xlsx+'/fyrment.xlsx',sheet_name='reactions',skiprows=[0]).fillna('').set_index('Rxn name')
		self.mets=pd.read_excel(path_fyrment_xlsx+'/fyrment.xlsx',sheet_name='metabolites').fillna('').set_index('Metabolite name')
		self.compartments=pd.read_excel(path_fyrment_xlsx+'/fyrment.xlsx',sheet_name='compartments').set_index('compartment_id').fillna('')
		print('\tLoad FOGs')
		self.fogs=self.get_fogs()
	def clean_fog_text(self,text_fog,keep_optional=False):
		delims='|&()'
		if keep_optional==False:
			delims+='?'
		for x in list(delims):
			text_fog=text_fog.replace(x,' ')
		return sorted(list(filter(None,set(text_fog.split(' ')))))
	def get_fyrment_fogs(self):
		text_fog=' '.join(self.rxns.FOG[1:])
		return self.clean_fog_text(text_fog)
	def get_reactions_from_fog(self,fog):
		return self.rxns[self.rxns.FOG.str.contains(fog)]
	def get_fogs_from_rxnid(self,rxnid,keep_optional):
		text_fog=self.rxns['FOG'][rxnid]
		return self.clean_fog_text(text_fog,keep_optional)



path_fyrment_xlsx='./fyrment'

fyrment=FYRMENT(path_fyrment_xlsx)



fyrment_fogs=[fog for fog in fyrment.fogs if fog in aybrah.df.index]

kla_fogs=[fog for fog in fyrment_fogs if len(aybrah.df['kla'][fog])>0]


