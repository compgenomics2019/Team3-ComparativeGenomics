#add shebang line
import argparse
import os
from os import listdir
from os.path import isfile, exists
import multiprocessing
import subprocess
import pandas as pf
import re
import sys

filenames=[]
filenames1=[]
filenames2=[]
input_dir=[]
dir_ref=[]

def check_files(path_to_file):
	#sys.exit if empty or not in correct format
	
	print("okay")
	#check is the file is in correct format
	#check if the file is not empty
	#check if the input directory exists
	return True

def parallelisation(i,j):
	#print(filenames)
	#command to run mummer ani
	c= "dnadiff -p tmp_ANI/"+filenames[i]+"to"+filenames[j] + " " + input_dir[0] + "/" +filenames[i] + "/scaffolds.fasta"+ " " + input_dir[0]+ "/" + filenames[j] + "/scaffolds.fasta" 
	check1=subprocess.run(c.split(), check=True, stdout=subprocess.PIPE, universal_newlines=True)
	#val=0
	with open("tmp_ANI/"+filenames[i]+"to"+filenames[j]+".report", "r") as f2:
		count=0
		count1=0
		for out in f2.readlines():
			out=out.strip("\n").strip(" ")
			if out[0:11]=="AvgIdentity" and count==0:
				for i in out.split(" "):
					if i!="AvgIdentity" and len(i)>1 and count1==0:
						val=i
						count1+=1
						return val
def parallelisation2(i,j):
	#command to run mummer ani

	c= "dnadiff -p tmp_ref/"+filenames1[i]+"to"+filenames2[j] + " " + input_dir[0] + "/" +filenames1[i] + "/scaffolds.fasta" + " " + dir_ref[0] + "/" + filenames2[j] 
	check1=subprocess.run(c.split(), check=True, stdout=subprocess.PIPE, universal_newlines=True)
	
	with open("tmp_ref/"+filenames1[i]+"to"+filenames2[j]+".report", "r") as f2:
		count=0
		count1=0
		for out in f2.readlines():
			out=out.strip("\n").strip(" ")
			if out[0:11]=="AvgIdentity" and count==0:
				for i in out.split(" "):
					if i!="AvgIdentity" and len(i)>1 and count1==0:
						val=i
						count1+=1
						return val

def compute_ani(filenames,thread,v,out):
	if v:
		print("Calculating Average Nucleotide Identity")
	lst_2=[]
	for i in range(0,len(filenames)):
		for j in range(i,len(filenames)):
			lst_2.append((i,j))

	exists = isfile(out)
	if exists:
		bash_cmd=["rm", out]
		subprocess.run(bash_cmd)

	if not os.path.exists("tmp_ANI"):
		command2= "mkdir tmp_ANI"
		subprocess.run(command2.split())

	results=[]
	p=multiprocessing.Pool(thread)
	results=p.starmap(parallelisation,lst_2)
	df= pf.DataFrame(0,columns= filenames, index= filenames)
	for i in range(0,len(results)):
		ind_1=lst_2[i][0]
		ind_2=lst_2[i][1]
		df.iloc[ind_1,ind_2]=results[i]
		df.iloc[ind_2,ind_1]=results[i]
	if v:
		print("ANI computation Completed")
		print("Writing to the output file")	

	if len(results)>0:
		writer=pf.ExcelWriter(out)
		df.to_excel(writer, sheet_name="Sheet1")
		writer.save()


def strain_identification(filenames2,thread,v,out_file,filenames1):
	if v:	
		print("Computation for strain Identification started")
	
	lst_3=[]
	for i in range(0,len(filenames1)):
		for j in range(0,len(filenames2)):
			lst_3.append((i,j))
	exists = isfile(out_file)
	if exists:
		bash_cmd=["rm", out_file]
		subprocess.run(bash_cmd)

	if not os.path.exists("tmp_ref"):
		command2= "mkdir tmp_ref"
		subprocess.run(command2.split())
	results=[]
	p=multiprocessing.Pool(thread)
	results=p.starmap(parallelisation2,lst_3)
	df= pf.DataFrame(0,columns= filenames2, index= filenames1)
	if len(results)>0:
		for i in range(0,len(results)):
			ind_1=lst_3[i][0]
			ind_2=lst_3[i][1]
			df.iloc[ind_1,ind_2]=results[i]
	if v:
		print("ANI computation between reference and isolates is Completed")
		print("Writing to the output file")	
	
	writer=pf.ExcelWriter(out_file)
	df.to_excel(writer, sheet_name="Sheet1")
	writer.save()

def compute_MLST(m,v, threshold): #compute_MLST(filenames,args.processes,args.verbose,args.outputANI)
	#Create wg schema:
	if v:	
		print("Creating folder for whole genome schema")
	"""if not os.path.exists("wgSchema"):
		c1= "mkdir wgSchema"
		subprocess.run(c1.split())
	if not os.path.exists("wg_alleleCall"):
		c2= "mkdir wg_alleleCall"
		subprocess.run(c2.split())"""
	
	command1="chewBBACA.py CreateSchema -i " + "tmp_genomes " + "-o" + " wgSchema" + " --cpu " + str(m) + " --ptf " + "Listeria_monocytogenes.trn"  
	subprocess.run(command1.split())
	#(get this file from comparative_genomics/nicole/prodigal_training/)

	#Allele calling for wg MLST:*
	"""command3= "chewBBACA.py AlleleCall -i " + "tmp_genomes" +  " -g" + " wgSchema " + "-o " + "wg_alleleCall " + "--cpu " + str(m) + " --ptf " + " Listeria_monocytogenes.trn" 

	#*Create cgMLST schema:*
	command4="chewBBACA.py ExtractCgMLST -i " + "/wg_alleleCall/results_alleles.tsv -o " + "cgSchema " + "-r " + "wg_alleleCall/RepeatedLoci.txt -p " + str(threshold) #0.95 (ask for user input for threshold?)
	
	subprocess.run(command1.split())
	subprocess.run(command3.split())
	subprocess.run(command4.split())
	if not os.path.exists("cgSchema/schema"):
		c1= "mkdir cgSchema/schema"
		subprocess.run(c1.split())
	
	#cd name/of/cgSchema/folder/
	#mkdir schema
	for files in os.listdir("wgSchema"):
		if files[-5:]=="fasta":
			cmd_1="cp " + files + " cgSchema/schema/"
	#cd schema
	#cat ../cgMLSTschema.txt | xargs -I {} cp name/of/wgSchema/folder/wgSchema/{} .*Allele calling for cg MLST:*
	command5="cd schema"
	command6="chewBBACA.py AlleleCall -i " + "tmp_genomes" + " -g " + "../cgSchema " + "-o " + " ../wg_alleleCall " + "--cpu " + str(m) + " --ptf Listeria_monocytogenes.trn"
	subprocess.run(command5.split())
	subprocess.run(command6.split())	"""

#def SNP():
#	command1="MakeFasta " + "tmp_SNP/out.txt " +  

def main():
	#define argparse arguments
	parser=argparse.ArgumentParser(description="Calculating the Genomic Distances", prog='pipeline.py', usage='%(prog)s [options]')
	parser.add_argument("-o", "--outputANI", required=False, type=str, nargs="?", const='output_ani.xlsx', default='output_ani.xlsx', help="provide output file name for ANI")
	parser.add_argument("-oMLST", required=False, type=str, help="provide output file name for cgMLST")
	parser.add_argument("-oSNP", required=False, type=str, help="provide output file name for SNP")
	parser.add_argument("-V", "--visualization", required=False, help="visualization required")
	parser.add_argument("-d", "--dir", required=False, type=str, help="path to input directory with the genome assemblies")
	parser.add_argument("-A", "--ANI", required=False, action="store_true",  help="compute ANI")
	parser.add_argument("-M", "--MLST", required=False, help="compute cgMLST", action="store_true")
	#parser.add_argument("-wM", "--wgMLST", required=False, help="compute cgMLST", action="store_true")
	parser.add_argument("-S", "--SNP", required=False, help="compute SNP analysis", action="store_true")
	parser.add_argument("-m", "--processes", required=False, type=int, default=1, help="provide the number of processes to launch")
	parser.add_argument("-v","--verbose", required=False, help="specify to know what is happening", action="store_true")
	#get names for strains	
	parser.add_argument("-R", "--strain_identification", required=False, default=False, help="perform strain Identification",action="store_true")
	parser.add_argument("-I", "--ref_dir", required=False, type=str, help="path to input directory with reference genomes")
	parser.add_argument("-Os", "--strain_out", required=False, type=str, nargs="?", default='output_ref.xlsx',const='output_ref.xlsx', help="provide output file name for comparison with reference genomes using ANI")
	parser.add_argument("-strains", "--strain_names", required=False, type=str, nargs="*", default='', help="name of the representative isolates")
	parser.add_argument("-t", "--threshold", required=False, type=int, default=0.95, help="provide the threshold for running cgMLST")
	#parser.add_argument("-f", "--force", required=False, help="specify if you want to overwrite existing files", action="store_true")
	
	args=parser.parse_args()

	#Input directory where genome assemblies are stored
	input_dir.append(args.dir)
	if not args.ANI or not args.SNP or args.MLST:
		print("Please provide the analysis that needs to be performed")
		print("Please check: pipeline.py --help")
		#sys.exit()
	for entry in args.strain_names:
		filenames1.append(entry)
	if len(input_dir[0])==0:
		print("Please provide directory with genome assemblies")
		sys.exit()

	if len(listdir(input_dir[0]))==0:
        	print("The input directory with genome assemblies is empty")
        	sys.exit()

	if args.strain_identification:
		dir_ref.append(args.ref_dir)
		if len(dir_ref[0])==0:
			print("Please provide directory with reference genomes")
			sys.exit()
		else:
			if len(listdir(args.ref_dir))==0:
        			print("The directory with reference genomes is empty")
        			sys.exit()
			if len(filenames1)==0:
				print("please perform pairwise ANI computation and provide representative genome ids from each cluster")
				sys.exit()
			else:
				#call commands to copy geome saccafolds for representatives to the all_reference folder
				lst_3=os.listdir(args.ref_dir)
				#check if all files are in correct format and non-empty
				for line2 in lst_3:
					with open(args.ref_dir+"/"+line2, "r") as f3:
						y=check_files(f3)
						if y:
							filenames2.append(line2)
				
				strain_identification(filenames2,args.processes,args.verbose,args.strain_out,filenames1) 	
							
	if args.ANI:
		lst_2=os.listdir(input_dir[0])
		for line in lst_2:
			with open(input_dir[0] +"/" + line +"/scaffolds.fasta", "r") as f2:
				x=check_files(f2)
				#if line not in filenames and len(f2.readlines())>0:
				if x:
					if line not in filenames:
						filenames.append(line)
		print("hi")
		compute_ani(filenames,args.processes,args.verbose,args.outputANI)
	if args.MLST or args.SNP:
		for folders in os.listdir(input_dir[0]):
			if not os.path.exists("tmp_genomes"):
				command2= "mkdir tmp_genomes"
				subprocess.run(command2.split())
			cmd1="cp -r " + input_dir[0]+"/"+folders + "/scaffolds.fasta " + "tmp_genomes/"
			cmd2="mv " + "tmp_genomes/scaffolds.fasta " + "tmp_genomes/" + folders + ".fasta"		
			subprocess.run(cmd1.split())
			subprocess.run(cmd2.split())

	if args.MLST:
		print("MLST")
		compute_MLST(args.processes,args.verbose,args.threshold)

	#if args.SNP:
	#	print("SNP") 

if __name__=="__main__":
    main()
