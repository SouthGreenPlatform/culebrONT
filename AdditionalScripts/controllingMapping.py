#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
# @package controllingMapping.py
# @author Francois Sabot

"""
	Controlling mapping data for a de novo assembly
	==========================
	:author:  François Sabot (thanks to Julie Orjuela Scripts help)
	:contact: francois.sabot@ird.fr
	:date: 04/04/2019
	:version: 0.1
	Script description
	------------------
	controllingMapping.py will take a BAM file issued from the mapping of Illumina data on a de novo assembly and try to catch the misassembled contigs based on local depth and local mapping errors
	-------
	>>> controllingMapping.py -i mapping.bam -o outputPrefix [-s windowSize -t tempLocation]
	
	Help Programm
	-------------
	information arguments:
		- \-h, --help
						show this help message and exit
		- \-v, --version
						display controllingMapping.py version number and exit
	Input mandatory infos for running:
		- \-i <filename>, --input <filename>
						BAM file issued from mapping (must be indexed)
		- \-o <filename>, --out <filename>
						Prefix of output files
		- \-s <windowSize>, --size <windowSize>
						window scan size in bases (Optional, default 1000)
		- \-t <tempLocation>, --size <tempLocation>
						Location for temp file, default /tmp
"""

##################################################
## Modules
##################################################
import sys, os, subprocess, re
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"

## Python modules
import argparse
from time import localtime, strftime
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##################################################
## Variables Globales
##################################################
version="0.1"
VERSION_DATE='04/04/2018'
debug="False"
#debug="True"

##################################################
## Functions
##################################################
def checkParameters (arg_list):
	# Check input related options
	if (not arg_list.inputFile):
		print ('Error: No input file defined via option -i/--input !' + "\n")
		parser.print_help()
		exit()
	if (not arg_list.outputPrefix):
		print ('Error: No output prefix file defined via option -o/--output !' + "\n")
		parser.print_help()
		exit()
		
def relativeToAbsolutePath(relative): 
	from subprocess import check_output
	if relative[0] != "/":			# The relative path is a relative path, ie do not starts with /
		command = "readlink -m "+relative
		absolutePath = subprocess.check_output(command, shell=True).decode("utf-8").rstrip()
		return absolutePath
	else:							# Relative is in fact an absolute path, send a warning
		absolutePath = relative;
		return absolutePath
	
def multicov(bedFile,bamFile,optionsList=("")):
	'''Commentaires'''
	optionsTxt = "".join(optionsList)
	command = "bedtools multicov " + optionsTxt + " -bams " + bamFile + " -bed " + bedFile
	#print (command)
	try:
		optionsTxt = str(OptionsTxt)
		resultCom = os.popen(command,"r")
	except NameError:
		resultCom = os.popen(command,"r")
	
	hashCov = {}
	for line in resultCom:
		#print(line)
		mainLine = line.strip()
		fields = mainLine.split("\t")
		start = int(fields[1])
		stop = int(fields[2])
		position = start + int(stop - start/2)
		try:
			hashCov[fields[0]].append((position,int(fields[3])))
		except KeyError:
			hashCov[fields[0]] = [(position,int(fields[3]))]
		#break
		#sys.exit()
	return hashCov	

##################################################
## Main code
##################################################
if __name__ == "__main__":
	# Parameters recovery
	parser = argparse.ArgumentParser(prog='controllingMapping.py', description='''This Program identifies the misassemblies in contigs using BAM, and uses BedTools2''')
	parser.add_argument('-v', '--version', action='version', version='You are using %(prog)s version: ' + version, help=\
						'display parsingSniffles.py version number and exit')
	filesreq = parser.add_argument_group('Input mandatory infos for running')
	filesreq.add_argument('-i', '--input', metavar="<filename>", required=True, dest = 'inputFile', help = 'vcf file')
	filesreq.add_argument('-o', '--out', metavar="<filename>", required=True, dest = 'outputPrefix', help = 'Prefix for the output files')
	
	parser.add_argument('-s', '--size', metavar="<sizeInBP>", required=False, default=1000, type=int, dest = 'windowSize', help = 'Window size for scanning')
	parser.add_argument('-t', '--temp', metavar="<diskLocation>", required=False, default="/tmp", dest = 'tempLocation', help = 'Location of temp directory')


	
	
	# Check parameters
	args = parser.parse_args()
	checkParameters(args)
	
	#Welcome message
	print("#################################################################")
	print("# 	Welcome in controllingMapping (Version " + version + ")	   	#")
	print("#################################################################")

	#Window size for scanning
	windowSize=args.windowSize
	#Temp location
	tempLocation=args.tempLocation


	# From relative to absolute paths 
	inputFile = relativeToAbsolutePath(args.inputFile)
	outputPrefix = relativeToAbsolutePath(args.outputPrefix)
	
	#open output handle
	outputFile = outputPrefix + "_position.csv"
	outputHandle = open(outputFile, "w")
	
	#Temp bed opening
	tempBed = tempLocation + "/temp.bed"
	tempHandle = open(tempBed, "w")

	#PySam import
	bamFile = pysam.AlignmentFile(inputFile,"rb");
	
	#Picking up infos on sequences
	sizes=bamFile.lengths
	names=bamFile.references
	seqInfo = {a:b for a,b in zip(names,sizes)}
	bamFile.close()
	
	for seq in seqInfo.keys():
		for window in range(1,seqInfo[seq],windowSize):
			start = str(window)
			stop = window + windowSize
			if stop > seqInfo[seq]:
				stop = seqInfo[seq] 
			stop = str(stop)
			tempHandle.write(seq + "\t" + start + "\t" + stop + "\n")
	
	depthTotal = multicov(tempBed, inputFile)
	options = ("-p")
	depthOk = multicov(tempBed, inputFile, options)
	
	#Creating the dataframe
	for contigs in depthTotal.keys():
		depthInfo={a:b for a,b in depthTotal[contigs]}
		okInfo={a:b for a,b in depthOk[contigs]}
		globalInfo={}
		for position in depthInfo.keys():
			globalInfo[position]=[depthInfo[position],okInfo[position]]
		df =pd.DataFrame.from_dict(globalInfo, orient='index',columns=['Depth','Ok'])
		df = df.sort_index()
		#print(df)
		somme = df['Depth'].sum()
		if somme > 0:
			plt.figure()
			ax = plt.gca()
			df.plot(kind='line',y='Depth', ax=ax,title=contigs)
			df.plot(kind='line',y='Ok',color='red',ax=ax)
			filename = contigs + ".png"
			plt.savefig(filename)
		else:
			continue
		
	sys.exit()
