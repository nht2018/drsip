'''
usage 

(bash)
python move_result.py --src_path mittelmann_ali/logs \
                          --dest_path MITTELMANN_ALI/logs \
                          --filename grb_ipm.log \
                          --copy true \
                          --rename grb_ipm.log

(python) 
from move_result import move_result

for method in ["grb_ipm",  "highs_ipm",  "pdlp", "grb_pspx", "grb_dspx"]:
    move_result("mittelmann_ali/results",
            "MITTELMANN_ALI/results",
            method + ".txt",
            copy=True)
    
for method in ["grb_ipm",  "highs_ipm", "grb_pspx", "grb_dspx"]:
    move_result("mittelmann_ali/logs",
            "MITTELMANN_ALI/logs",
            method + ".log",
            copy=True,
            rename=method + ".log")

'''



import os
import shutil
import argparse

def move_result(src_dir, 
                dest_dir,
                filename="grb_ipm.txt", 
                instance_list=None, 
                copy=True, #if true, then copy the files; otherwise move the files
                rename=""):#give the files new names
    
    if instance_list is None:
        instance_list = os.listdir(src_dir)
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)
    for instance in instance_list:
        src_path = os.path.join(src_dir, instance, filename)
        if rename != "":
            dest_path = os.path.join(dest_dir, instance, rename)
        else:
            dest_path = os.path.join(dest_dir, instance, filename)
        if not os.path.exists(src_path):
            print(src_path, " not exist")
            continue
        if not os.path.exists(os.path.join(dest_dir, instance)):
            os.mkdir(os.path.join(dest_dir, instance))
        if copy:
            shutil.copyfile(src_path, dest_path)
        else:
            shutil.move(src_path, dest_path)

def delete_result(src_dir, 
                filename="grb_ipm.txt", 
                instance_list=None):
    
    if instance_list is None:
        instance_list = os.listdir(src_dir)
    for instance in instance_list:
        src_path = os.path.join(src_dir, instance, filename)
        if not os.path.exists(src_path):
            print(src_path, " not exist")
            continue
        os.remove(src_path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--src_dir", type=str, default="MITTELMANN_ALI/result")
    parser.add_argument("--dest_dir", type=str, default="MITTELMANN_ALI/results")
    parser.add_argument("--filename", type=str, default="grb_dspx.txt")
    parser.add_argument("--copy", type=lambda x: x.lower() == 'true', default=True)
    parser.add_argument("--rename", type=str, default="")
    args = parser.parse_args()

    move_result(args.src_dir, args.dest_dir, 
        filename=args.filename, 
        copy=args.copy,
        rename=args.rename)


if __name__=="__main__":
    main()