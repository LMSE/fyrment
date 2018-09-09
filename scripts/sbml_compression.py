
import pandas as pd
import gzip
import os


"""
files were written in xml locally.
In future could compress immediatly to GitHub folder
As the total storage of all the GENRE's is approaching 1 GB

"""


taxonomy=pd.read_csv('taxonomy_level_oids.txt',sep='\t')
#taxonomy.drop(46)
drop_these=[index for oid in ['cpr','ani'] for index in taxonomy[taxonomy['oids']==oid].index.tolist()]
taxonomy=taxonomy.drop(drop_these)

taxon_order=zip(*sorted([(taxonomy.level.tolist().index(taxon),taxon) for taxon in set(taxonomy.level.tolist())]))[1]




for index,row in taxonomy.iterrows():
	level=taxonomy['level'][index]
	level_order=taxon_order.index(level)
	name=taxonomy['name'][index]
	readable=taxonomy['readable'][index]
	print readable
	for extension in ['xml','xlsx']:
		print '\t'+extension
		path_genre_dir = '/Users/kcorreia/Dropbox/PhD/04_manuscripts/MOPUS03_AYbRAHAM/'+extension+'/'+str(level_order).zfill(2)+'_'+level.lower()
		path_genre_file=path_genre_dir+'/i'+readable+'.'+extension
		path_genre_gzip_dir='./genre/'+extension+'/'+path_genre_file.split('/')[-2]
		if not os.path.exists(path_genre_gzip_dir):
			os.makedirs(path_genre_gzip_dir)
		f_in = open(path_genre_file)
		f_out = gzip.open(path_genre_gzip_dir+'/i'+readable+'.'+extension+'.gz', 'wb')
		f_out.writelines(f_in)
		f_out.close()
		f_in.close()

