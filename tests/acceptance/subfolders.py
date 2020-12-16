from __future__ import print_function
from os import listdir
import os

#refactor code to intergrate with test_accept_chunk.py
root_dir="/home/xuyanting/results"
new_dir="/home/xuyanting/results_new"
files_per_folder=1000

file_names=[]
for input_file in listdir(root_dir):
    #make a new sub dir for every 1000 files
    file_names.append(input_file)
print(len(file_names))

folders = []
for fold_id in range(0,int(len(file_names)/files_per_folder)):
   folders.append(new_dir+"/"+str(fold_id))
   os.mkdir(new_dir+"/"+str(fold_id))
counter=0

print("number of fodlers",len(folders))
for fold in folders:
   print("folder", fold)
   for index in range(counter,counter+files_per_folder):
       print(root_dir+"/"+file_names[index])
       print(fold+"/"+file_names[index])
       os.rename(root_dir+"/"+file_names[index], fold+"/"+file_names[index])
   counter += files_per_folder
