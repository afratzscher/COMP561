'''
FILE: 1000GPfetch.py
PURPOSE: connects to 1000GP ftp and fetches SNP data and creates txt file
INPUT: none
OUTPUT: error code (-1 if error)
'''

import config
import os
import pandas as pd
from pathlib import Path
import subprocess

# def autosomes(filepath):
# 	config.__FILENAME__ = "1000G_chr" + (str(int(config.__CHR__)) + ".vcf")
# 	return


def run_commands(*commands):
	os.system(' ; '.join(commands))
def combine():
	raw = read_vcf.read_vcf(config.__FILEPATH__+'raw_chr' + str(config.__CHR__) + "_" + str(config.__START__) + "-" 
			+ str(config.__END__) + ".vcf")
	chrom = str(config.__CHR__)

	dbSNP = pd.read_csv(config.__FILEPATH__+'dbSNP_chr'+ chrom + "_" + str(config.__START__) + "-" 
		+ str(config.__END__) + ".vcf", sep='\t', header=None) # ISSUE HERE
	dbSNP.columns = ['POS', 'ID', 'REF', 'ALT']

	for i in range(0, len(raw.index)):
		# if NOT in dbSNP, dont replace (remove later)
		if not (dbSNP.loc[dbSNP['POS'] == raw.loc[i]['POS']]).empty: # check if pos match
			row = dbSNP.loc[dbSNP['POS'] == raw.loc[i]['POS']].values[0]
			if ((raw.iloc[i,3] == row[2]) and (raw.iloc[i,4] == row[3])): # only replace if REF/ALT match
				raw.iloc[i,2] = row[1] # ID
				raw.iloc[i,3] = row[2] # REF
				raw.iloc[i,4] = row[3] # ALT
	df = raw[raw.ID != '.'] # remove if not in dbSNP
	df.to_csv((config.__FILEPATH__ + config.__FILENAME__), sep="\t", mode='a', index=False)

# for P3 grch38 b154
def makeCommandsP3_38_154(name, ftp, cmds, xy):
	if cmds == "":
		cmds = []
	#define ftp
	cmds.append(ftp+name)

	#define dbSNP ftp -> GRCH38 BUILD 154 
	baseUrl = 'dbSNPFTP=ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/'
	refVersion = 'GCF_000001405.38.gz'
	cmds.append(baseUrl+refVersion)

	#command to get data file
	baseName =  config.__FILENAME__[len('1000G_'):] 
	cmds.append('bcftools view "$ftp" -I -r ' + str(config.__CHR__) + ":" 
			+ str(config.__START__) + "-" + str(config.__END__) 
			+ " > " + config.__FILEPATH__ + "raw_" + baseName)
	cmds.append('bcftools query "$dbSNPFTP" -f "%POS\t%ID\t%REF\t%ALT\n" -r ' + config.__CHRVERSION__
					+ ":" + str(config.__START__) + "-" + str(config.__END__) + "> " + config.__FILEPATH__ + "dbSNP_" + baseName)
	
	# clean up
	cmds.append('rm ' + name + '.tbi')
	cmds.append('rm ' + refVersion + '.tbi')
	return cmds
def getData(filepath):
	createFolder(filepath) # create folder if doesnt exist

	for i in range(1, 22):
		rawName = "raw_" + config.__FILENAME__[len('1000G_'):] 
		if not Path(config.__FILEPATH__ + rawName).is_file():
			ftp = "ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/"
			vcfgzName = "ALL.chr" + str(i) + "_GRCh38.genotypes.20170504.vcf.gz"
			cmd = makeCommandsP3_38_154(vcfgzName, ftp, "", "")
			run_commands(*cmd)
		
		#then combine raw and dbSNP
		# combine()
		break
		# autosomes(filepath)

def main():
	print("*****STARTING SELECTION*****")
	config.__FOLDERPATH__ = os.getcwd()[:-(len('src/'))]
	errCode = getData(config.__FOLDERPATH__)
	
if __name__ == '__main__':
	main()