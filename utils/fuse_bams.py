import os
import glob
import argparse

parser = argparse.ArgumentParser("Generate a list of all the files in a directory")

parser.add_argument("-d", "--directory", type=str, required=True, help="The directory to search")

args = parser.parse_args()
dir_of_folders = args.directory

items = os.listdir(dir_of_folders)
folders = [item for item in items if os.path.isdir(os.path.join(dir_of_folders, item))]

for folder in folders:
    print(os.path.join(dir_of_folders, folder))
