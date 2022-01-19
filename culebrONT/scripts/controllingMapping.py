#!/usr/bin/env python3
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

"""

##################################################
## Modules
##################################################
import sys, os, subprocess, re
current_dir = os.path.dirname(os.path.abspath(__file__))+"/"

## Python modules
import click
from pathlib import Path
from datetime import datetime
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


@click.command("controlling_mapping", short_help='This Program identifies the misassemblies in contigs using BAM, and uses BedTools2',
			   context_settings={'help_option_names': ('-h', '--help'),"max_content_width":800}, no_args_is_help=True)
@click.version_option(version, '--version', '-v')
@click.option('--input', '-i', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True), required=True, help='BAM file issued from mapping (must be indexed)')
@click.option('--out', '-o', type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True, resolve_path=True), required=True, help='Prefix for the output files')
@click.option('--size', '-s', type=int, required=False, default=1000, show_default=True, help='window scan size in bases')
@click.option('--temp', '-t', type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True), required=False, default="/tmp", show_default=True, help='Location of temp directory')
def main(input, out, size, temp):
	"""
	\b
	This Program identifies the misassemblies in contigs using BAM, and uses BedTools2
	controllingMapping.py will take a BAM file issued from the mapping of Illumina data on a de novo assembly
	and try to catch the misassembled contigs
	based on local depth and local mapping errors
	"""
	inputFile = input
	outputPrefix = out
	windowSize = size
	tempLocation = temp
	start_time = datetime.now()

	# Welcome message
	click.secho(f"""{"#" * 80}\n#{Path(sys.argv[0]).stem + " " + version:^78}#\n{"#" * 80}\nStart time: {start_time:%d-%m-%Y at %H:%M:%S}\nCommande line run: {" ".join(sys.argv)}\n""", fg="green")

	# open output handle
	outputFile = Path(outputPrefix).joinpath("_position.csv")
	outputHandle = open(outputFile, "w")

	# PySam import
	bamFile = pysam.AlignmentFile(inputFile, "r");

	# Picking up infos on sequences
	sizes = bamFile.lengths
	names = bamFile.references
	seqInfo = {a: b for a, b in zip(names, sizes)}
	bamFile.close()

	# Temp bed opening
	with open(tempLocation + "/temp.bed", "w") as tempHandle:
		for seq in seqInfo.keys():
			for window in range(1, seqInfo[seq], windowSize):
				start = str(window)
				stop = window + windowSize
				if stop > seqInfo[seq]:
					stop = seqInfo[seq]
				stop = str(stop)
				tempHandle.write(seq + "\t" + start + "\t" + stop + "\n")

	depthTotal = multicov(tempBed, inputFile)
	options = ("-p")
	depthOk = multicov(tempBed, inputFile, options)

	# Creating the dataframe
	for contigs in depthTotal.keys():
		depthInfo = {a: b for a, b in depthTotal[contigs]}
		okInfo = {a: b for a, b in depthOk[contigs]}
		globalInfo = {}
		for position in depthInfo.keys():
			globalInfo[position] = [depthInfo[position], okInfo[position]]
		df = pd.DataFrame.from_dict(globalInfo, orient='index', columns=['Depth', 'Ok'])
		df = df.sort_index()
		# print(df)
		somme = df['Depth'].sum()
		if somme > 0:
			plt.figure()
			ax = plt.gca()
			df.plot(kind='line', y='Depth', ax=ax, title=contigs)
			df.plot(kind='line', y='Ok', color='red', ax=ax)
			filename = contigs + ".png"
			plt.savefig(filename)
		else:
			continue

	sys.exit()
##################################################
## Main code
##################################################
if __name__ == "__main__":
	main()
