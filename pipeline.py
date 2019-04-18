#!/usr/bin/env python
import argparse
import os
from os import listdir
from os.path import isfile, exists
import multiprocessing
import subprocess
import pandas as pf
import re
import sys
from os.path import expanduser

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
						return 100-val
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
						return 100-val

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

def compute_MLST(m,v, threshold,output_cgMLST): #compute_MLST(filenames,args.processes,args.verbose,args.outputANI)
	
	if v:	
		print("Running whole genome and core genome MLST")
	if not os.path.exists("tmp_MLST"):
		comand1= "mkdir tmp_MLST"
		subprocess.run(comand1.split())
		comd_tmp="mv Listeria_monocytogenes.trn tmp_MLST/"
		subprocess.run(comd_tmp.split())
	
	command1="chewBBACA.py CreateSchema -i " + "tmp_genomes " + "-o" + " tmp_MLST/wgSchema" + " --cpu " + str(m) + " --ptf " + " tmp_MLST/Listeria_monocytogenes.trn"  

	#Allele calling for wg MLST:*
	command3= "chewBBACA.py AlleleCall -i " + "tmp_genomes" +  " -g" + " tmp_MLST/wgSchema " + "-o " + "tmp_MLST/wg_alleleCall " + "--cpu " + str(m) + " --ptf " + " tmp_MLST/Listeria_monocytogenes.trn" 
	subprocess.run(command1.split())
	subprocess.run(command3.split())
	tmp_file=os.listdir("tmp_MLST/wg_alleleCall")[0]

	command4="chewBBACA.py ExtractCgMLST -i " + "tmp_MLST/wg_alleleCall/" + tmp_file + "/results_alleles.tsv -o " + "tmp_MLST/cgSchema " + "-r " + "tmp_MLST/wg_alleleCall/" + tmp_file + "/RepeatedLoci.txt -p " + str(threshold) #0.95 (ask for user input for threshold?)
	
	subprocess.run(command4.split())
	if not os.path.exists("tmp_MLST/cgSchema/schema"):
		c1= "mkdir tmp_MLST/cgSchema/schema"
		subprocess.run(c1.split())
	with open("tmp_MLST/cgSchema/cgMLSTschema.txt", "r") as f1:
		for line in f1.readlines():
			line=line.strip("\n").strip(" ")
		#for files in os.listdir("tmp_MLST/wgSchema"):
			cmd_1="cp " + "tmp_MLST/wgSchema/"+line + " tmp_MLST/cgSchema/schema"
			subprocess.run(cmd_1.split())
	cmd_ee="cp -r " + "tmp_MLST/wgSchema/short" + " tmp_MLST/cgSchema/schema/"
	subprocess.run(cmd_ee.split())
	#command5="cd tmp_MLST/cgSchema/schema"
	#subprocess.run(command5.split())
	#command_ee="pwd ."
	#subprocess.run(command_ee.split())
	#command6="chewBBACA.py AlleleCall -i " + "tmp_MLST/cgSchema/schema/../../../tmp_genomes" + " -g " + "tmp_MLST/cgSchema/schema/../../cgSchema " + "-o " + " tmp_MLST/cgSchema/schema/../../wg_alleleCall " + "--cpu " + str(m) + " --ptf tmp_MLST/cgSchema/schema/../../Listeria_monocytogenes.trn"
	#print(command6.split())
	#current=os.getcwd()
	command6="chewBBACA.py AlleleCall -i " + "tmp_genomes" + " -g " + "tmp_MLST/cgSchema/schema " + "-o " + "tmp_MLST/cg_alleleCall " + "--cpu " + str(m) + " --ptf " + "tmp_MLST/Listeria_monocytogenes.trn"
	print(command6.split())
	#subprocess.check_call(command6.split(), cwd="tmp_MLST/cgSchema/schema")
	subprocess.run(command6.split())
	tmp_file2=os.listdir("tmp_MLST/cg_alleleCall")[0]
	command_term= "mv " + "tmp_MLST/cg_alleleCall/" + tmp_file2 + "/results_alleles.tsv" + " " + output_cgMLST
	subprocess.run(command_term.split())

def compute_SNP(v,outdir, enviro):
	if v:	
		print("Running KSNP")
	if not os.path.exists("tmp_SNP"):
		comand1= "mkdir tmp_SNP"
		subprocess.run(comand1.split())
		comand2="mkdir outdir"
		subprocess.run(comand2.split())
	#command1="MakeFasta " + "tmp_SNP/out.txt " +  "tmp_SNP/fastainput"
 	#/path/to/input_list.txt /path/to/fastainput
	#Kchooser /path/to/fastainput
	current=os.getcwd()
	for file in os.listdir("tmp_genomes"):
		with open("tmp_SNP/temp.txt", "a") as f2:
			f2.write(current+"tmp_genomes/"+file+"\n")
	home = expanduser("~")
	#command2= home+"/miniconda2/envs/"+enviro+"/bin"+"/kSNP3 -in " + current+"/tmp_SNP/temp.txt" + " -outdir " + current +"/"+ outdir + " -k 19 -vcf -ML -NJ" 
	#subprocess.run(command2.split())
	command2="kSNP3 -in " + current+"/tmp_SNP/temp.txt" + " -outdir " + current +"/"+ outdir + " -k 19 -vcf -ML -NJ"
	print(command2)	
	subprocess.check_call(command2.split(), cwd=home+"/miniconda2/envs/"+enviro+"/bin")
	
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
	parser.add_argument("-outsnp", "--out_SNP", required=False, type=str, nargs="?", const='output_SNP.xlsx', default='output_SNP.xlsx', help="provide output file name for SNP analysis")
	parser.add_argument("-outmlst", "--out_MLST", required=False, type=str, nargs="?", const='output_MLST.xlsx', default='output_MLST.xlsx', help="provide output file name for cgMLST analysis")
	parser.add_argument("-e", "--environment", required=False, type=str, default="", help="provide the name of the conda environment")
	#parser.add_argument("-f", "--force", required=False, help="specify if you want to overwrite existing files", action="store_true")
	
	args=parser.parse_args()

	#Input directory where genome assemblies are stored
	input_dir.append(args.dir)
	if not args.ANI or not args.SNP or not args.MLST:
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
		compute_MLST(args.processes,args.verbose,args.threshold. args.out_MLST)

	if args.SNP: 
		if len(args.environment)>0:
			compute_SNP(args.verbose,args.out_SNP,args.environment)
		else:
			print("Please provide name of the conda enviornment")
			sys.exit()

if __name__=="__main__":
    main()
