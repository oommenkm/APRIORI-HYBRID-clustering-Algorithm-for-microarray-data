"""
This Module will perform the T-Test on the data values in the file given

Writen on: 20th Feb 2012

Author: Oommen Mathew
"""

import math
import sys

def prepareValue(x):
    x = re.sub(r'"|,|positive|negative')
    parts = x.split('e+')
    return float(parts[0]) * math.exp(int(parts[1]))
    
class TTest:
    
    def __init__(self, file_name):
        self.lns = None
        self.pos_cont=0
        self.neg_cont=0
        FH = open(file_name,'r')
        self.lns = FH.readlines()
        FH.close()
        
        self.t_test_value = []
        self.pos_value_sum = []
        self.neg_value_sum = []
        self.pos_values = []
        self.neg_values = []
        
        for i in range(0,1999):
            self.pos_values.append([])
            self.neg_values.append([])
            self.pos_value_sum.append(0)
            self.neg_value_sum.append(0)
            
    def calculateTTestValues(self):
        for ln in self.lns:
            flds=ln.split('\t')
            if(flds[1999].find('positive')>=0):
                self.pos_cont=self.pos_cont+1
                for i in range(1999):
                    cur_fld_val = prepareValue(flds[i])
                    self.pos_values[i].append(cur_fld_val)
                    self.pos_value_sum[i] = self.pos_value_sum[i] + cur_fld_val
                    
            else:
                self.neg_cont=self.neg_cont+1
                for i in range(1999):
                    cur_fld_val = prepareValue(flds[i])
                    self.neg_values[i].append(cur_fld_val)
                    self.neg_value_sum[i]=self.neg_value_sum[i]  + cur_fld_val
        
        positive_mean = []
        negative_mean = []
        for i in range(1999):
            positive_mean.append(self.pos_value_sum[i]/self.pos_cont)
            negative_mean.append(self.neg_value_sum[i]/self.neg_cont)
            
        mean_dif_value = []
        for i in range(1999):
            mean_dif_value.append(positive_mean[i]-negative_mean[i])
            
        pos_standard_dev = []
        neg_standard_dev = []
        for i in range(1999):
            pos_square_dif_sum = 0
            for j in range(self.pos_cont):
                x_mean = self.pos_values[i][j] - positive_mean[i]
                
                pos_square_dif_sum = pos_square_dif_sum + (x_mean * x_mean)
            
            pos_standard_dev.append(math.sqrt((pos_square_dif_sum / (self.pos_cont - 1))))
            
            neg_square_dif_sum = 0
            for j in range(self.neg_cont):
                x_mean = self.neg_values[i][j] - negative_mean[i]
                neg_square_dif_sum = neg_square_dif_sum + (x_mean * x_mean)
            
            neg_standard_dev.append(math.sqrt((neg_square_dif_sum / (self.neg_cont - 1))))
        
        sum_std_dev = []
        sqrt_std_dev = []
        for i in range(1999):
            sum_std_dev = pos_standard_dev[i] + neg_standard_dev[i]
            sqrt_std_dev = math.sqrt(sum_std_dev)
            self.t_test_value.append(float(mean_dif_value[i]/sqrt_std_dev))
            #print i,t_test_value

    def get_tops(self, num):
        tops = []
        top_ids = []
            
        for i in xrange(num):
            tops.append(-999)
            top_ids.append(None)
        
        for i in range(1999):
            for j in range(num):
                if(float((self.t_test_value[i]) >= tops[j])):
                    for p in range(num-1, j, -1):
                        tops[p] = tops[p-1]
                        top_ids[p] = top_ids[p - 1]
                    tops[j] = float(self.t_test_value[i])
                    top_ids[j] = i
                    break
        if -999 in top_ids:
            final_top = top_ids[0:top_ids.index(-999)]
        else:
            final_top = top_ids
        return final_top
    
    def print_tops(self):
        for i in xrange(num):
            print tops[i], "-", top_ids[i]
            
def main():
    t_obj = TTest(sys.argv[1])
    t_obj.calculateTTestValues()
    top_id_list = t_obj.get_tops(sys.argv[2])
    print top_id_list
    
if __name__ == '__main__':
    main()