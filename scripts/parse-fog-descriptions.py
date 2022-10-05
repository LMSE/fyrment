

df=pd.read_excel('aybrah.xlsx').set_index('FOG').fillna('')


columns_description=['HOG', 'Protein ID', 'sce Locus', 'Features','Protein description', 'Parent', 'Neofunc, Subfunct', 'References', 'Source', 'Suggested Analysis','SGD Description', 'PomBase Description', 'AspGD Description']
readable=pd.DataFrame(columns=columns_description)

for fog,row in df.iterrows():
	print(fog)
	row_extracted=[row[x] for x in columns_description]
	if len(list(filter(None,row_extracted[1:8]+row_extracted[9:])))==0:
		continue
	readable.loc[fog]=row_extracted


readable.to_csv('./aybrah/descriptions_fog.txt',sep='\t')
