#!/usr/bin/env python

import glob
import gc
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.linecharts import HorizontalLineChart

final_result_dict = {}

for infile in glob.glob('seven_Whole*.csv' ):
    print "current file is: " + infile
    cur_no_genes = infile.split("_Whole_")[1].split(".csv")[0]
    
    print "Analyzing ", cur_no_genes, "...Genes"

    whole_file_name = "seven_Whole_" + cur_no_genes + ".csv"
    pos_file_name = "seven_Positive_" + cur_no_genes + ".csv"
    neg_file_name = "seven_Negative_" + cur_no_genes + ".csv"
    
    wholeFH = open(whole_file_name, 'r')
    posFH = open(pos_file_name, 'r')
    negFH = open(neg_file_name, 'r')
    
    #whole_file_lines = wholeFH.readlines()
    #pos_file_lines = posFH.readlines()
    #neg_file_lines = negFH.readlines()
    
    pos_rules = []
    for ln in posFH.readlines():
        pos_rules.append(ln.split(",")[0])
    
    print "Len of pos_rules : ", len(pos_rules)
    
    neg_rules = []
    for ln in negFH.readlines():
        neg_rules.append(ln.split(",")[0])
    
    print "Len of neg_rules : ", len(neg_rules)
    
    pos_count_in_whole = 0
    neg_count_in_whole = 0
    whole_file_lines = 0
    for ln in wholeFH.readlines():
        whole_file_lines = whole_file_lines + 1
        cur_rule = ln.split(",")[0]
        if cur_rule in pos_rules:
            pos_count_in_whole = pos_count_in_whole + 1
        if cur_rule in neg_rules:
            neg_count_in_whole = neg_count_in_whole + 1

    cur_res = {}
    cur_res["pos_count"] = pos_count_in_whole
    cur_res["neg_count"] = neg_count_in_whole
    cur_res["pos_per"] = round(pos_count_in_whole/(float(pos_count_in_whole + neg_count_in_whole)) * 100.00, 2)
    cur_res["neg_per"] = round(neg_count_in_whole/(float(pos_count_in_whole + neg_count_in_whole)) * 100.00, 2)
    uncat = whole_file_lines - (pos_count_in_whole + neg_count_in_whole)
    cur_res["uncat"] = uncat
    cur_res["uncat_per"] = round(uncat/float(whole_file_lines) * 100.00, 2)
    
    final_result_dict[cur_no_genes] = cur_res
    
    finalOutPutFH = open("finalPutput.txt", 'w')
    finalOutPutFH.write(str(final_result_dict))
    finalOutPutFH.close()
    
    wholeFH.close()
    posFH.close()
    negFH.close()
    del(cur_no_genes)
    del(pos_rules)
    del(neg_rules)
    del(cur_res)
    
    gc.collect()
    
normal_pos = 22.00 / 62.00 * 100.00
normal_neg = 40.00 / 62.00 * 100.00

#print "\n\n+---------|----------|--------|",
#for ky in final_result_dict.keys():
#    print "--------|",
#print "--------+"
#print "|Type     |Normal %|Count     |",
#for ky in final_result_dict.keys():
#    print ky, "   |",
#print "\n+---------|----------|--------|",
#for ky in final_result_dict.keys():
#    print "--------|",
#print "--------+"

print "No. of Genes  | Pos count  | Neg Count | Pos %   | Neg %   | Uncategorized Count | Uncategorized Per |"

for cnt in final_result_dict.keys():
    print cnt,
    print final_result_dict[cnt]["pos_count"], "|",
    print final_result_dict[cnt]["neg_count"], "|",
    print final_result_dict[cnt]["pos_per"], "|",
    print final_result_dict[cnt]["neg_per"], "|",
    print final_result_dict[cnt]["neg_per"], "|",
    print final_result_dict[cnt]["uncat"], "|",
    print final_result_dict[cnt]["uncat_per"], "|"

#print "+---------|----------|--------|--------+"

