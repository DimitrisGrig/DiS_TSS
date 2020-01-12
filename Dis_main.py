import os
import sys
from scipy.stats import kurtosis, skew
from numpy import trapz, array
from statistics import mean
import pickle
import pandas as pd
from scipy.signal import find_peaks
import subprocess

def read_file(file):
	with open(file) as f:
		current = []
		for line in f:  # read rest of lines
			current.append([x for x in line.strip().split("\t")])
	return(current)

def split_file_name(filename, delim):
	tmp = filename.split(delim)  # Split filename
	return(tmp[0])

def _args_(input_file):
	_args = input_file.split("/")
	_args.pop(0)
	path = "/"+'/'.join(map(str,_args[0:len(_args)-1]))
	filename = _args[len(_args)-1]
	name = split_file_name(filename, ".")
	return (path, filename, name)

def create_dir(path):
	if not os.path.isdir(path):
		os.system("mkdir " + path)

def alter_read_1bp(filename, outpufile_pos, outpufile_neg):
	with open(filename, 'r') as infile:
		with open(outpufile_pos, 'w') as outfile_pos:
			with open(outpufile_neg, 'w') as outfile_neg:
				for line in infile:
					line = line.rstrip().split("\t")
					if line[5] == "+":
						# For '+' strand
						line_pos = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(line[0], line[1], str(int(line[1])+1), line[3], line[4], line[5])
						outfile_pos.write(line_pos)
					elif line[5] == "-":
						# For '-' strand
						line_neg = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(line[0], str(int(line[2])-1), line[2], line[3], line[4], line[5])
						outfile_neg.write(line_neg)

def split_strand(file, outpufile_pos, outpufile_neg):
	with open(file, 'r') as infile:
		with open(outpufile_pos, 'w') as outfile_pos:
			with open(outpufile_neg, 'w') as outfile_neg:
				for line in infile:
					line = line.rstrip().split("\t")
					if line[5] == "+":
						# For '+' strand
						line_pos = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(line[0], line[1], line[2], line[3], line[4], line[5])
						outfile_pos.write(line_pos)
					elif line[5] == "-":
						# For '-' strand
						line_neg = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(line[0], line[1], line[2], line[3], line[4], line[5])
						outfile_neg.write(line_neg)

def lookup_table(filename):
	lookup = {}
	depth = []
	with open(filename) as f:
		for idx, line in enumerate(f):
			e0, e1 = line.split()[0], line.split()[1]
			lookup[str(e0)+','+str(e1)] = idx   # use e0,e1 as key, and idx as value
			depth.append(line.rstrip().split()[2])
	return lookup, depth

def predict_svm(features_file, work_dir, features_list):
	# Importing the dataset
	dataset = pd.read_csv(features_file, header=None)

	X = dataset.iloc[:, features_list].values
	X = X.astype(float)

	print("Feature Scaling...")

	# load scaler
	script_dir = os.path.dirname(os.path.abspath(__file__))
	sc = pickle.load(open(script_dir+"/model/scaler.pkl", 'rb'))

	# Feature Scaling
	X_test = sc.transform(X)

	# Load SVM model
	classifier = pickle.load(open(script_dir+"/model/model.mdl", 'rb'))

	# Predicting results
	y_pred_proba = classifier.predict_proba(X_test)[:,1]

	return(y_pred_proba)

def clustering(input_file, read_quality, number_of_threads, map_distance, tmp_folder,tpm_threshold):
	[work_dir, filename, name] = _args_(input_file)
	out = work_dir+"/"+name+"_quality.bam"
	print("Clustering...")
	os.system("samtools view -@ "+str(number_of_threads)+" -q "+str(read_quality)+" -b "+input_file+" > "+out)
	current = out

	# Count total reads
	command = "samtools view -c "+current+" -@ "+str(number_of_threads)
	result = subprocess.check_output(command, shell=True).rstrip()
	total_reads = 1000000 / int(result.decode("utf-8"))

	out = work_dir+"/"+name+"_quality.bed"
	os.system("bedtools bamtobed -i "+current+" > "+out)
	current = out

	out = work_dir+"/"+name+"_quality_1bp.bed"
	with open(current, 'r') as infile:
		with open(out, 'w') as outfile:
			for line in infile:
				line = line.split("\t")
				# Change 3rd column of file ==> 2nd column ++1
				if line[5].rstrip() == "+":
					#For '+' strand
					new_line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(line[0], line[1], int(line[1])+1, line[3], line[4], line[5])
				elif line[5].rstrip() == '-':
					#For '-' strand
					new_line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(line[0], int(line[2])-1, line[2], line[3], line[4], line[5])

				#Print in File
				outfile.write(new_line)
	current = out

	out = work_dir+"/"+name+"_quality_1bp_sorted.bed"
	os.system("sort -k1,1 -k2,2n -V -s "+current+" > "+out)
	current = out

	final_clusters = work_dir+"/"+name+"_clusters.bed"
	os.system("bedtools merge -i "+current+" -s -d "+str(map_distance)+" -c 1,6 -o count,distinct > "+final_clusters)

	out = work_dir+"/"+name+"_clusters_tmp.bed"
	os.system("awk  -F \"\t\" -v OFS=\"\t\" '{ if ( $4 != \"1\" ) { print $1,$2,$3,$(NF+1) = \"tc_tss_\"(NR-1) 1 ,"+ str(total_reads)+"*$4, $5; } }' "+final_clusters+" > "+out)
	current = out

	out = work_dir+"/"+name+"_clusters_final.bed"
	os.system("awk  -F \"\t\" -v OFS=\"\t\" '{ if ( $5 > "+str(tpm_threshold)+" ) { print $0; } }' "+current+" > "+out)
	return(out)
#


############################ MAIN ############################
############################ INPUTS ############################
# Script dir
script_dir = os.path.dirname(os.path.abspath(__file__))

args = list(sys.argv)
if len(sys.argv) > 0 and len(sys.argv) < 13:
	if '-i' in args: 
		input_file = args[args.index('-i')+1] 
	else:
		print("Input file not specified.\n")
		print("Examle usage of tool:\n"
	       +"python Dis_main.py -i <input_file> [OPTIONS]\n\n"
	       +"-i         - CAGE-seq, bam file, (provide full path e.g. /home/user/data.bam)\n"
	       +"-clusters  - Bed file containing clusters. If set to \"0\" the tool will auto estimate the clusters. Default is 0 (auto).\n"
	       +"-@         - Number of theads to use. Default is 1.\n"
	       +"-MAPQ      - Skip alignments with map quality lower than <-MAPQ>. Default is 10.\n"
	       +"-MAPDist   - Minimum distance for CTSS to be merged. Default is 25.\n"
	       +"-tpm       - Clusters expression below <-tpm> are expelled. Default is 1.\n")
		sys.exit()
	if '-clusters' in args: 
		clusters_path = args[args.index('-clusters')+1] 
	else: 
		clusters_path = 0
	if '-@' in args: 
		number_of_threads = args[args.index('-@')+1] 
	else:
		 number_of_threads = 1
	if '-MAPQ' in args: 
		read_quality = args[args.index('-MAPQ')+1]
	else:
		read_quality = 10
	if '-MAPDist' in args:
		map_distance = args[args.index('-MAPDist')+1]
	else:
		map_distance = 25
	if '-tpm' in args:
		tpm_threshold = args[args.index('-tpm')+1]
	else:
		tpm_threshold = 1
else:
	print("DiS-TSS takes exactly 1 to 6 arguments.\n")
	print("Examle usage of tool:\n"
	       +"python Dis_main.py -i <input_file> [OPTIONS]\n\n"
	       +"-i         - CAGE-seq, bam file (provide full path e.g. /home/user/data.bam)\n"
	       +"-clusters  - Bed file containing clusters. If set to \"0\" the tool will auto estimate the clusters. Default is 0 (auto).\n"
	       +"-@         - Number of theads to use. Default is 1.\n"
	       +"-MAPQ      - Skip alignments with map quality lower than <-MAPQ>. Default is 10.\n"
	       +"-MAPDist   - Minimum distance for CTSS to be merged. Default is 25.\n"
	       +"-tpm       - Clusters expression below <-tpm> are expelled. Default is 1.\n")

# Determine working directory and file name
[work_dir, filename, name] = _args_(input_file)

tmp_folder = work_dir+"/tmp"
create_dir(tmp_folder)

# Utilize clustering
if clusters_path == 0 or clusters_path == "0":
	clusters_path = clustering(input_file, read_quality, number_of_threads, map_distance, tmp_folder, tpm_threshold)

model = script_dir+"/model/model.mdl"
features_list = [9,11,8,10,6,7]

chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7',
			'chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15',
			'chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

# Sort clusters file
clusters_path_name = os.path.splitext(clusters_path)[0]
os.system("sort -k1,1 -V -s "+clusters_path+" > "+clusters_path_name+"_sorted.bed")

# Read clusters file
clusters = read_file(clusters_path_name+"_sorted.bed")

# Create export folder
create_dir(work_dir+"/files")

print("Sorting bam...")
# Compute Bam coverage #
if not os.path.isfile(work_dir+"/files/"+name+"_sorted.bam"):
	os.system("samtools sort -o "+work_dir+"/files/"+name+"_sorted.bam -T "+tmp_folder+" -O bam -@ "+str(number_of_threads)+" "+work_dir+"/"+filename)
	filename = name+"_sorted.bam"
	os.system("samtools index "+work_dir+"/files/"+filename+" -@ "+str(number_of_threads))
else:
	filename = name+"_sorted.bam"

print("BamToBed...")
# Bam To Bed
if not os.path.isfile(work_dir+"/files/"+name+".bed"):
	os.system("bedtools bamtobed -i "+work_dir+"/files/"+filename+" > "+work_dir+"/files/"+name+".bed")

print("Changing tag size...")
# Change tag size to first nucleotide 
if not os.path.isfile(work_dir+"/files/"+name+"_1bp_pos.bed") or not os.path.isfile(work_dir+"/files/"+name+"_1bp_neg.bed"):
	alter_read_1bp(work_dir+"/files/"+name+".bed", work_dir+"/files/"+name+"_1bp_pos.bed", work_dir+"/files/"+name+"_1bp_neg.bed")

# Bed To Bam
if not os.path.isfile(work_dir+"/files/"+name+"_1bp_pos.bam") or not os.path.isfile(work_dir+"/files/"+name+"_1bp_neg.bam"):
	os.system("bedtools bedtobam -ubam -i "+work_dir+"/files/"+name+"_1bp_pos.bed -g /mnt/raid0/dgrigor/annotation/hg38.fa.fai > "+work_dir+"/files/"+name+"_1bp_pos.bam")
	os.system("bedtools bedtobam -ubam -i "+work_dir+"/files/"+name+"_1bp_neg.bed -g /mnt/raid0/dgrigor/annotation/hg38.fa.fai > "+work_dir+"/files/"+name+"_1bp_neg.bam")

	# Sort
	os.system("samtools sort -o "+work_dir+"/files/"+name+"_1bp_pos_sorted.bam -T "+tmp_folder+" -O bam -@ "+str(number_of_threads)+" "+work_dir+"/files/"+name+"_1bp_pos.bam")
	os.system("samtools sort -o "+work_dir+"/files/"+name+"_1bp_neg_sorted.bam -T "+tmp_folder+" -O bam -@ "+str(number_of_threads)+" "+work_dir+"/files/"+name+"_1bp_neg.bam")

print("Splitting tags on strand...")
# Split clusters_path file on strand
if not os.path.isdir(work_dir+"/clusters"):
	create_dir(work_dir+"/clusters")
	split_strand(clusters_path, work_dir+"/clusters/positions_pos.bed", work_dir+"/clusters/positions_neg.bed")

print("Computing coverage...", end="")
# Compute coverage
if not os.path.isdir(work_dir+"/coverage"):
	create_dir(work_dir+"/coverage")
	os.system("samtools depth -a -Q "+str(read_quality)+" -b "+work_dir+"/clusters/positions_pos.bed "+work_dir+"/files/"+name+"_1bp_pos_sorted.bam > "+work_dir+"/coverage/"+name+"_1bp_pos_cov.bebgraph")
	os.system("samtools depth -a -Q "+str(read_quality)+" -b "+work_dir+"/clusters/positions_neg.bed "+work_dir+"/files/"+name+"_1bp_neg_sorted.bam > "+work_dir+"/coverage/"+name+"_1bp_neg_cov.bebgraph")
print("DONE")

lookup_pos, depth_pos = lookup_table(work_dir+"/coverage/"+name+"_1bp_pos_cov.bebgraph")
lookup_neg, depth_neg = lookup_table(work_dir+"/coverage/"+name+"_1bp_neg_cov.bebgraph")

print("Calculating features")
# Open file where features of clusters are exported/stored
with open(work_dir+"/"+name+"_clusters_features.bed", 'w') as f:
	with open(work_dir+"/"+name+"_clusters_features.csv", 'w') as f_csv:
		# For each cluster (for each line in clusters file)
		for cluster in clusters:
			if cluster[0] in chr_list:
				if cluster[5] == "-":
					search_key = str(cluster[0]) + ',' + str(int(cluster[1])+1)
					start = lookup_neg[search_key]

					search_key = str(cluster[0]) + ',' + str(cluster[2])
					end = lookup_neg[search_key]

					signal = depth_neg[start:end+1]
					strand = -1
				elif cluster[5] == "+":
					search_key = str(cluster[0]) + ',' + str(int(cluster[1])+1)
					start = lookup_pos[search_key]

					search_key = str(cluster[0]) + ',' + str(cluster[2])
					end = lookup_pos[search_key]

					signal = depth_pos[start:end+1]
					strand = 1

				# cast "signal" list to int
				cur_signal = list(map(int, signal))

				# Compute features
				_peak_length = len(cur_signal)
				_peak_hight = max(cur_signal)   # Max depth
				_kurt = kurtosis(cur_signal)
				_skew = skew(cur_signal)
				peaks, properties = find_peaks(cur_signal, prominence=0, height=0, width=0)
				if len(peaks) > 1:
					_mean = mean(peaks)
				else:
					_mean = 0

				print(*cluster[0:6], _peak_length, _peak_hight, round(_kurt,2), round(_skew,2), round(_mean,2), len(peaks), strand, sep = "\t", file=f)
				print(*cluster[0:6], _peak_length, _peak_hight, round(_kurt,2), round(_skew,2), round(_mean,2), len(peaks), strand, sep = ",", file=f_csv)

print("Predicting...")
proba = predict_svm(work_dir+"/"+name+"_clusters_features.csv", work_dir, features_list)

clusters_ft = read_file(work_dir+"/"+name+"_clusters_features.bed")

with open(work_dir+"/"+name+"_prob_clusters.bed", 'w') as f1:
	for idx, cluster in enumerate(clusters_ft):
		print(*cluster[0:3], str(cluster[3])+"||"+str(cluster[4]), round(proba[idx], 2), cluster[5], sep = "\t", file=f1)

if os.path.isdir(tmp_folder):
	os.system("rm -d "+tmp_folder)

# Remove middle files
if os.path.isfile(work_dir+"/H9_quality.bed"):
	os.system("rm "+work_dir+"/H9_quality.bed")

if os.path.isfile(work_dir+"/H9_quality_1bp.bed"):
	os.system("rm "+work_dir+"/H9_quality_1bp.bed")

if os.path.isfile(work_dir+"/H9_quality_1bp_sorted.bed"):
	os.system("rm "+work_dir+"/H9_quality_1bp_sorted.bed")   

if os.path.isfile(work_dir+"/H9_quality.bam"):
	os.system("rm "+work_dir+"/H9_quality.bam")

if os.path.isfile(work_dir+"/H9_clusters_tmp.bed"):
	os.system("rm "+work_dir+"/H9_clusters_tmp.bed")

if os.path.isfile(work_dir+"/H9_clusters_final_sorted.bed"):
	os.system("rm "+work_dir+"/H9_clusters_final_sorted.bed")

if os.path.isfile(work_dir+"/H9_clusters.bed"):
	os.system("rm "+work_dir+"/H9_clusters.bed")

if os.path.isfile(work_dir+"/H9_clusters_features.csv"):	
	os.system("rm "+work_dir+"/H9_clusters_features.csv")

if os.path.isfile(work_dir+"/H9_clusters_features.bed"):
	os.system("rm "+work_dir+"/H9_clusters_features.bed")

if os.path.isdir(work_dir+"/files"):
	os.system("rm -r "+work_dir+"/files")

if os.path.isdir(work_dir+"/clusters"):
	os.system("rm -r "+work_dir+"/clusters")

print("Finished")
#