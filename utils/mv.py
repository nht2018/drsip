'''
Author: Hantao Nie (nht@pku.edu.cn)
Date: 2023-10-28 19:06:42
LastEditors: Hantao Nie (nht@pku.edu.cn)
LastEditTime: 2023-10-28 19:32:12
Description: 

Copyright (c) 2023, Hantao Nie, Peking University. 
'''
from move_result import move_result, delete_result

for method in ["drspf", "SDPT3", "mosek", "sedumi", "gurobi"]:
    move_result(src_dir="/Users/niehantao/Desktop/software/ssn_sqlp/results/CBLIB/1e-6",
                dest_dir="/Users/niehantao/Desktop/software/conic_programming_framework/results/CBLIB",
                filename=method + ".txt", copy=False)
