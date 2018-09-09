[  if len(aybrah[oid][fog[:8]])>0 ]


rxnid='BIOMASS_PROTEIN_M_iMM904_TRNA_1G'

for fog in complex_essential:
	print fog,aybrah[oid][fog]



for fname in $(cat models.txt) ; do
  echo $fname
  curl -F file=@"$fname" -F output=xml -F offcheck=u http://sbml.org/validator/ > "$fname""_validation.xml"
done


fname=Strain_Hanseniaspora_valbyensis_NRRL_Y-1626_AYbRAH_0.1.8_GSMM.xml
curl -F file=@"$fname" -F output=xml -F offcheck=u http://sbml.org/validator/ > "$fname""_validation.xml"


