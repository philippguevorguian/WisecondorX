from tqdm import tqdm
import pysam
import os
import glob
import argparse

parser = argparse.ArgumentParser("Generate a list of all the files in a directory")

parser.add_argument("-d", "--directory", type=str, required=True, help="The directory to search")

parser.add_argument("-o", "--output", type=str, required=True, help="The directory to output to")

args = parser.parse_args()
dir_of_folders = args.directory
output_dir = args.output

items = os.listdir(dir_of_folders)
folders = [item for item in items if os.path.isdir(os.path.join(dir_of_folders, item))]

for folder in tqdm(folders):
    folder_path = os.path.join(dir_of_folders, folder)
    tel_bam = glob.glob(folder_path + '/**/tel.align.bam', recursive=True)[0]
    base_bam = glob.glob(folder_path + '/**/base.align.bam',recursive = True)[0]
    print(folder)
    folder_name = str(folder)
    bam_name = folder_name + ".bam"
    merged_bam_file = os.path.join(output_dir,bam_name)
    os.system(f"samtools merge {merged_bam_file} {base_bam} {tel_bam}")
    pysam.index(merged_bam_file)
    
