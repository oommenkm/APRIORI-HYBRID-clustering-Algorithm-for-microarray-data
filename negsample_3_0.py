
import sys
import math
import gc
import TTest

def prepareValue(x):
    x = re.sub(r'"|,|positive|negative')
    parts = x.split('e+')
    return float(parts[0]) * math.exp(int(parts[1]))

def prune_genes(whole_support, indiv_gene_value_support):
    pruned_dis = []
    for c_k in indiv_gene_value_support.keys():
        if indiv_gene_value_support[c_k] > whole_support:
            pruned_dis.append(c_k)
            
    return pruned_dis

(i, j, k, l, m, n) = (-1, -1, -1, -1, -1, -1)

if sys.argv[1]:
    top_limit = int(sys.argv[1])
else:
    top_limit = 500

try:
    t_obj = TTest.TTest('colon data.txt')
    t_obj.calculateTTestValues()
    top_id_list = t_obj.get_tops(top_limit)
    
    print "Selected Gene Ids : ", top_id_list[0:10], "..."
    
    sel_gene_count = len(top_id_list)
    
    FH = open('colon data.txt','r')
    #FH = open('C:\Users\Merry\Desktop\colon data.txt','r')
    value_cont = 0
    lns = FH.readlines()

    FH.close()

    sum = 0

    row_value = {}
    
    for ln in lns:
        flds=ln.split('\t')
        row_col_value = []
        if(flds[1999].find('negative')>=0):
            for i in top_id_list:
                if not row_value.has_key(i):
                    row_value[i] = []
                value_cont = value_cont + 1
                cur_value = prepareValue(flds[i])
                row_col_value.append(cur_value)
                sum = sum + cur_value
                row_value[i].append(cur_value)
            
    for ky in row_value.keys():
        if ky in top_id_list[0:10]:
            print ky, '---', row_value[ky][0:3], "..."
    
    del(row_col_value)
    del(i)
    del(ln)
    del(flds)
    del(lns)
    
    gc.collect()
        
    value_mean = sum / value_cont
    
    value_dis = {}
    
    count_of_value_ones = 0
    for cur_key in row_value.keys():
        p_val = row_value[cur_key]
        value_dis_row = []
        for i in range(len(p_val)):
            cur_value = p_val[i]
            if(cur_value >= value_mean):
                cur_value = 1
                count_of_value_ones = count_of_value_ones + 1
                value_dis_row.append(cur_value)        
            else:
                cur_value = 0
                value_dis_row.append(cur_value)
                
        value_dis[cur_key] = value_dis_row
        
    whole_support = count_of_value_ones / (value_cont * 1.00)
    
    for c_k in row_value.keys():
        if c_k in top_id_list[0:10]:
            print c_k, value_dis[c_k][0:10]
    
    indiv_gene_value_support = {}    
    for cur_k in row_value.keys():
        indiv_gene_value_support[cur_k] = value_dis[cur_k].count(1) / 40.00
        
    for c_k in row_value.keys():
        if c_k in top_id_list[0:10]:
            print c_k, indiv_gene_value_support[c_k]
    
    pruned_value_dis = prune_genes(whole_support, indiv_gene_value_support)
    print len(pruned_value_dis)
    print pruned_value_dis[0:10]
    
    pruned_gene_list_sorted = pruned_value_dis
    pruned_gene_list_sorted.sort()
        
    gc.collect()
        
    trans1_conf_value_vals = []
    conf_value_vals  = []
    
    file_name_suffix = "_Negative_" + str(top_limit)
    
    twoFH = open("two" + file_name_suffix + ".csv", 'w')
    
    for i in range(len(pruned_gene_list_sorted) - 1):
        for j in range(i + 1, len(pruned_gene_list_sorted)):
            asso_count = 0
            indv_count = 0
            two_conf_list = []
            first_gene_id = pruned_gene_list_sorted[i]
            second_gene_id = pruned_gene_list_sorted[j]
            for cur_sample in range(40):
                if (value_dis[first_gene_id][cur_sample] == 1  and  value_dis[second_gene_id][cur_sample] == 1):
                    asso_count = asso_count + 1
                    two_conf_list.append('1')
                else:
                    two_conf_list.append('0')
                if (value_dis[first_gene_id][cur_sample] == 1):
                    indv_count = indv_count + 1
                    
            conf_value = asso_count/float(indv_count)
            if(conf_value > 0.5):
                comb_disc_val_str = ",".join(two_conf_list)
                twoFH.write(str(i) + "-" + str(j) + "," + comb_disc_val_str + "\n")
                twoFH.flush()
    
    twoFH.close()
    print "Finished combinaions of 2"
    
    print "Starting Combinations of 3"
    
    twoFH = open("two" + file_name_suffix + ".csv", 'r')    
    threeFH = open("three" + file_name_suffix + ".csv", 'w')
    
    curLine = ''
    prevLine = ''
    beforePrevLine = ''
    lineNumber = 0
    for ln in twoFH.readlines():
        ln = ln.replace('\n', '')
        ln = ln.replace('\r', '')
        lineNumber = lineNumber + 1
        beforePrevLine = prevLine
        prevLine = curLine
        curLine = ln
        if lineNumber < 3:
            continue
        #print '.',
        asso_count = 0
        indv_count = 0
        three_conf_list = []
        curLineFields = curLine.split(',')
        prevLineFields = prevLine.split(',')
        beforePrevLineFields = beforePrevLine.split(',')
        #for row_num in range(value_cont-1):
        for row_num in range(len(curLineFields) - 1):
            #print row_num, ',',
            cur_conf_val = '0'
            if(beforePrevLineFields[row_num + 1] == '1'  and  prevLineFields[row_num + 1] == '1' and curLineFields[row_num + 1] == '1'):
                asso_count = asso_count + 1
                cur_conf_val = '1'
            if (beforePrevLineFields[row_num + 1] == '1'  and  prevLineFields[row_num + 1] == '1'):
                indv_count = indv_count + 1
                
            three_conf_list.append(cur_conf_val)
            
        if indv_count * asso_count > 0:
            conf_value = asso_count/float(indv_count)
            if(conf_value > 0.5):
                comb_disc_val_str = ",".join(three_conf_list)
                threeFH.write(beforePrevLineFields[0] + "|" +
                              prevLineFields[0] + "|" +
                              curLineFields[0] + "," +
                              comb_disc_val_str + "\n")
                threeFH.flush()
    
    threeFH.close()
    twoFH.close()
    
    print "Finished combinaions of 3"
    
    print "Starting Combinations of 4"
    
    threeFH = open("three" + file_name_suffix + ".csv", 'r')    
    fourFH = open("four" + file_name_suffix + ".csv", 'w')
    
    curLine = ''
    prevLine = ''
    beforePrevLine = ''
    threeBeforePrevLine = ''
    lineNumber = 0
    for ln in threeFH.readlines():
        ln = ln.replace('\n', '')
        ln = ln.replace('\r', '')
        lineNumber = lineNumber + 1
        threeBeforePrevLine = beforePrevLine
        beforePrevLine = prevLine
        prevLine = curLine
        curLine = ln
        if lineNumber < 4:
            continue
        #print '.',
        asso_count = 0
        indv_count = 0
        four_conf_list = []
        curLineFields = curLine.split(',')
        prevLineFields = prevLine.split(',')
        beforePrevLineFields = beforePrevLine.split(',')
        threeBeforePrevLineFields = threeBeforePrevLine.split(',')
        #for row_num in range(value_cont-1):
        for row_num in range(len(curLineFields) - 1):
            #print row_num, ',',
            cur_conf_val = '0'
            if(threeBeforePrevLineFields[row_num + 1] == '1'  and  beforePrevLineFields[row_num + 1] == '1'  and  prevLineFields[row_num + 1] == '1' and curLineFields[row_num + 1] == '1'):
                asso_count = asso_count + 1
                cur_conf_val = '1'
            if (threeBeforePrevLineFields[row_num + 1] == '1'  and  beforePrevLineFields[row_num + 1] == '1'  and  prevLineFields[row_num + 1] == '1'):
                indv_count = indv_count + 1
            
            four_conf_list.append(cur_conf_val)
            
        if indv_count * asso_count > 0:
            conf_value = asso_count/float(indv_count)
            if(conf_value > 0.5):
                comb_disc_val_str = ",".join(four_conf_list)
                fourFH.write(threeBeforePrevLineFields[0] + ">" +
                             beforePrevLineFields[0] + ">" +
                             prevLineFields[0] + ">" +
                             curLineFields[0] + "," + 
                             comb_disc_val_str + "\n")
                fourFH.flush()
    
    fourFH.close()
    threeFH.close()
    
    print "Finished combinaions of 4"

    print "Starting Combinations of 5"
    
    fourFH = open("four" + file_name_suffix + ".csv", 'r')    
    fiveFH = open("five" + file_name_suffix + ".csv", 'w')
    
    curLine = ''
    prevLine = ''
    beforePrevLine = ''
    threeBeforePrevLine = ''
    fourBeforePrevLine = ''
    lineNumber = 0
    for ln in fourFH.readlines():
        ln = ln.replace('\n', '')
        ln = ln.replace('\r', '')
        lineNumber = lineNumber + 1
        fourBeforePrevLine =  threeBeforePrevLine
        threeBeforePrevLine = beforePrevLine
        beforePrevLine = prevLine
        prevLine = curLine
        curLine = ln
        if lineNumber < 5:
            continue
        #print '.',
        asso_count = 0
        indv_count = 0
        four_conf_list = []
        curLineFields = curLine.split(',')
        prevLineFields = prevLine.split(',')
        beforePrevLineFields = beforePrevLine.split(',')
        threeBeforePrevLineFields = threeBeforePrevLine.split(',')
        fourBeforePrevLineFields = fourBeforePrevLine.split(',')
        #for row_num in range(value_cont-1):
        for row_num in range(len(curLineFields) - 1):
            #print row_num, ',',
            cur_conf_val = '0'
            if(fourBeforePrevLineFields[row_num + 1] == '1'  and
               threeBeforePrevLineFields[row_num + 1] == '1'  and
               beforePrevLineFields[row_num + 1] == '1'  and
               prevLineFields[row_num + 1] == '1' and
               curLineFields[row_num + 1] == '1'):
                    asso_count = asso_count + 1
                    cur_conf_val = '1'
            if (fourBeforePrevLineFields[row_num + 1] == '1'  and
                threeBeforePrevLineFields[row_num + 1] == '1'  and
                beforePrevLineFields[row_num + 1] == '1'  and
                prevLineFields[row_num + 1] == '1'):
                    indv_count = indv_count + 1
                    
            four_conf_list.append(cur_conf_val)
            
        if indv_count * asso_count > 0:
            conf_value = asso_count / float(indv_count)
            if(conf_value > 0.5):
                comb_disc_val_str = ",".join(four_conf_list)
                fiveFH.write(fourBeforePrevLineFields[0] + "<" +
                             threeBeforePrevLineFields[0] + "<" +
                             beforePrevLineFields[0] + "<" +
                             prevLineFields[0] + "<" +
                             curLineFields[0] + "," + 
                             comb_disc_val_str + "\n")
                fiveFH.flush()
    
    fiveFH.close()
    fourFH.close()
    
    print "\nFinished combinaions of 5"
    
    print "Starting Combinations of 6"
    
    fiveFH = open("five" + file_name_suffix + ".csv", 'r')    
    sixFH = open("six" + file_name_suffix + ".csv", 'w')
    
    curLine = ''
    prevLine = ''
    beforePrevLine = ''
    threeBeforePrevLine = ''
    fourBeforePrevLine = ''
    fiveBeforePrevLine = ''
    lineNumber = 0
    for ln in fiveFH.readlines():
        ln = ln.replace('\n', '')
        ln = ln.replace('\r', '')
        lineNumber = lineNumber + 1
        fiveBeforePrevLine =  fourBeforePrevLine
        fourBeforePrevLine =  threeBeforePrevLine
        threeBeforePrevLine = beforePrevLine
        beforePrevLine = prevLine
        prevLine = curLine
        curLine = ln
        if lineNumber < 6:
            continue
        #print '.',
        asso_count = 0
        indv_count = 0
        four_conf_list = []
        curLineFields = curLine.split(',')
        prevLineFields = prevLine.split(',')
        beforePrevLineFields = beforePrevLine.split(',')
        threeBeforePrevLineFields = threeBeforePrevLine.split(',')
        fourBeforePrevLineFields = fourBeforePrevLine.split(',')
        fiveBeforePrevLineFields = fiveBeforePrevLine.split(',')
        #for row_num in range(value_cont-1):
        for row_num in range(len(curLineFields) - 1):
            #print row_num, ',',
            cur_conf_val = '0'
            if(fiveBeforePrevLineFields[row_num + 1] == '1'  and
               fourBeforePrevLineFields[row_num + 1] == '1'  and
               threeBeforePrevLineFields[row_num + 1] == '1'  and
               beforePrevLineFields[row_num + 1] == '1'  and
               prevLineFields[row_num + 1] == '1' and
               curLineFields[row_num + 1] == '1'):
                    asso_count = asso_count + 1
                    cur_conf_val = '1'
            if (fiveBeforePrevLineFields[row_num + 1] == '1'  and
                fourBeforePrevLineFields[row_num + 1] == '1'  and
                threeBeforePrevLineFields[row_num + 1] == '1'  and
                beforePrevLineFields[row_num + 1] == '1'  and
                prevLineFields[row_num + 1] == '1'):
                    indv_count = indv_count + 1
                    
            four_conf_list.append(cur_conf_val)
            
        if indv_count * asso_count > 0:
            conf_value = asso_count / float(indv_count)
            if(conf_value > 0.5):
                comb_disc_val_str = ",".join(four_conf_list)
                sixFH.write(fourBeforePrevLineFields[0] + "~" +
                             threeBeforePrevLineFields[0] + "~" +
                             beforePrevLineFields[0] + "~" +
                             prevLineFields[0] + "~" +
                             curLineFields[0] + "," + 
                             comb_disc_val_str + "\n")
                sixFH.flush()
    
    sixFH.close()
    fiveFH.close()
    
    print "\nFinished combinaions of 6"

    print "Starting Combinations of 7"
    
    sixFH = open("six" + file_name_suffix + ".csv", 'r')    
    sevenFH = open("seven" + file_name_suffix + ".csv", 'w')
    
    curLine = ''
    prevLine = ''
    beforePrevLine = ''
    threeBeforePrevLine = ''
    fourBeforePrevLine = ''
    fiveBeforePrevLine = ''
    sixBeforePrevLine = ''
    lineNumber = 0
    for ln in sixFH.readlines():
        ln = ln.replace('\n', '')
        ln = ln.replace('\r', '')
        lineNumber = lineNumber + 1
        sixBeforePrevLine =  fiveBeforePrevLine
        fiveBeforePrevLine =  fourBeforePrevLine
        fourBeforePrevLine =  threeBeforePrevLine
        threeBeforePrevLine = beforePrevLine
        beforePrevLine = prevLine
        prevLine = curLine
        curLine = ln
        if lineNumber < 7:
            continue
        #print '.',
        asso_count = 0
        indv_count = 0
        four_conf_list = []
        curLineFields = curLine.split(',')
        prevLineFields = prevLine.split(',')
        beforePrevLineFields = beforePrevLine.split(',')
        threeBeforePrevLineFields = threeBeforePrevLine.split(',')
        fourBeforePrevLineFields = fourBeforePrevLine.split(',')
        fiveBeforePrevLineFields = fiveBeforePrevLine.split(',')
        sixBeforePrevLineFields = sixBeforePrevLine.split(',')
        #for row_num in range(value_cont-1):
        for row_num in range(len(curLineFields) - 1):
            #print row_num, ',',
            cur_conf_val = '0'
            if(sixBeforePrevLineFields[row_num + 1] == '1'  and
               fiveBeforePrevLineFields[row_num + 1] == '1'  and
               fourBeforePrevLineFields[row_num + 1] == '1'  and
               threeBeforePrevLineFields[row_num + 1] == '1'  and
               beforePrevLineFields[row_num + 1] == '1'  and
               prevLineFields[row_num + 1] == '1' and
               curLineFields[row_num + 1] == '1'):
                    asso_count = asso_count + 1
                    cur_conf_val = '1'
            if (sixBeforePrevLineFields[row_num + 1] == '1'  and
                fiveBeforePrevLineFields[row_num + 1] == '1'  and
                fourBeforePrevLineFields[row_num + 1] == '1'  and
                threeBeforePrevLineFields[row_num + 1] == '1'  and
                beforePrevLineFields[row_num + 1] == '1'  and
                prevLineFields[row_num + 1] == '1'):
                    indv_count = indv_count + 1
                    
            four_conf_list.append(cur_conf_val)
            
        if indv_count * asso_count > 0:
            conf_value = asso_count / float(indv_count)
            if(conf_value > 0.5):
                comb_disc_val_str = ",".join(four_conf_list)
                sevenFH.write(fourBeforePrevLineFields[0] + "~" +
                             threeBeforePrevLineFields[0] + "~" +
                             beforePrevLineFields[0] + "~" +
                             prevLineFields[0] + "~" +
                             curLineFields[0] + "," + 
                             comb_disc_val_str + "\n")
                sevenFH.flush()
    
    sevenFH.close()
    sixFH.close()
    
    print "\nFinished combinaions of 7"
    
except Exception, e:
    print "Excepion Occurred:  ", e
    print i, j, k, l, m, n
    raise