## Title     : More detailed integration of information 
## Objective : Result of more detailed integration of information
## Created by: fpeng
## Created on: 2024/01/12


import sys
import pandas as pd
import re
import numpy as np


## Extract the required data from the original sort file and output the "blast_sortpy_result" file in a fixed format ##
def sites_info(name):
    site_results = ""
    with open(name,"r") as result:
        blast_all = result.readlines()
    site_results += "genome1_name" +"\t" + "genome1_site" + "\t"  + "genome2_name" +"\t" + "genome2_site" + "\t" + "gap or overlap" + "\t" +"reads_name" + "\n"
    for i in range(0,len(blast_all),2):
        genome1_name = blast_all[i].split("\t")[1]
        genome1_site = blast_all[i].split("\t")[9]
        genome2_name = blast_all[i+1].split("\t")[1]
        genome2_site = blast_all[i+1].split("\t")[8]
        gap = abs(int(blast_all[i].split("\t")[7])-int(blast_all[i+1].split("\t")[6]))
        reads_name = blast_all[i].split("\t")[0]
        
      	
        ## sort HPV and human so that human is in the first column ##
        if "HPV" in genome2_name:
            site_results += genome1_name +"\t" + genome1_site + "\t"  + genome2_name +"\t" + genome2_site + "\t" + str(gap) + "\t" + reads_name + "\n"
        else:
            site_results += genome2_name + "\t" + genome2_site + "\t" + genome1_name + "\t" + genome1_site + "\t" + str(
                gap) + "\t" + reads_name + "\n"
        #print(blast_all[i],end="")
  	
    ## write to the file ##
    with open("blast_sortpy_result","w") as f2:
        print(site_results,file=f2)
        f2.close()
    return site_results


## Barcode will be presented for each cell ##
def sites_cluster(info):
    cluster_pre = pd.read_csv(info,sep='\t')
    #print(info.split("\n")[1].split("\t")[0])
    #print(cluster_pre)
    
  	## Dataframe is sorted by genome name and locus ##
    cluster_sort = cluster_pre.sort_values(by=["genome1_name","genome1_site"],ascending=(True,True))
    cluster_sort = cluster_sort.reset_index(drop=True)
    #print(cluster_sort)
    cell = []
    for i in cluster_sort['reads_name']:
        barcode = re.findall(r"_(.+?)_",i)
        cell.append(barcode)
    cell = list(np.ravel(cell))
    cluster_sort['cell'] = cell

    cluster_sort_bycell = cluster_sort.sort_values(by=["cell"], ascending=(True))
    cluster_sort_bycell = cluster_sort_bycell.reset_index(drop=True)
    col = cluster_sort_bycell.pop('cell')
    cluster_sort_bycell.insert(loc=0,column='cell',value=col)
    #print(cluster_sort_bycell)
    
    ## View all the elements and display all the dataframe ##
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    cluster_sort_bycell_setindex = cluster_sort_bycell.sort_values(by=["cell","genome1_name","genome1_site",genome2_name","genome2_site"], ascending=(True,True,True,True,True))
    cluster_sort_bycell_setindex = cluster_sort_bycell_setindex.reset_index(drop=True)
  	##Get Rid of the gap no need #
    cluster_sort_bycell_setindex = cluster_sort_bycell_setindex.drop(['gap or overlap'],axis=1)
    cluster_sort_bycell_setindex.to_csv('cluster_sort_bycell_setindex',sep='\t',index=False,header=True,escapechar=' ')


## Statistical“sites_cluster_reslut” for the same barcode ##
def total_sum(sort_after):
    # siteswith open("sites_cluster_reslut.txt","w"):

    with open(sort_after,"r") as sort_txt:
        sorted_txt = sort_txt.readlines()
    #print(sorted_txt[1])
    sites_cluster_result = ''
    last_cotent = sorted_txt[1]
    num = 1
  	## The data of reads were calculated ##
    reads_contain = [last_cotent.split("\t")[5]]
    for i in sorted_txt[2:]:
        if i == sorted_txt[-1]:
            reads_contain.append(sorted_txt[-1].split("\t")[5])
            #print(last_cotent.split("\t")[0],"\t",last_cotent.split("\t")[1],"\t",last_cotent.split("\t")[2],"\t",last_cotent.split("\t")[3],"\t",last_cotent.split("\t")[4],"\t",num+1,"\t",reads_contain)
            sites_cluster_result += last_cotent.split("\t")[0]+"\t"+last_cotent.split("\t")[1]+"\t"+last_cotent.split("\t")[2]+"\t"+last_cotent.split("\t")[3]+"\t"+last_cotent.split("\t")[4]+"\t"+str(num+1)+"\t"+str(reads_contain)+"\n"
        elif i.split("\t")[0] == last_cotent.split("\t")[0] and i.split("\t")[1] == last_cotent.split("\t")[1] and abs(int(i.split("\t")[2])-int(last_cotent.split("\t")[2])) <=20 and i.split("\t")[3] == last_cotent.split("\t")[3] and abs(int(i.split("\t")[4])-int(last_cotent.split("\t")[4])) <=20:
            num += 1
            reads_contain.append(i.split("\t")[5])
            #print(i)
        else:
            sites_cluster_result += last_cotent.split("\t")[0] + "\t" + last_cotent.split("\t")[1] + "\t" + last_cotent.split("\t")[2] + "\t" + last_cotent.split("\t")[3] + "\t" + last_cotent.split("\t")[4] + "\t" + str(num) + "\t" + str(reads_contain) + "\n"
            #print(last_cotent.split("\t")[0],"\t",last_cotent.split("\t")[1],"\t",last_cotent.split("\t")[2],"\t",last_cotent.split("\t")[3],"\t",last_cotent.split("\t")[4],"\t",num,"\t",reads_contain)
            num = 1
            last_cotent = i
            reads_contain = [last_cotent.split("\t")[5]]
        #print(i)
  	## open the result file and write ##
    with open('sites_cluster_reslut','w') as f:
        f.write("cell_barcode"+"\t"+"human_genome"+"\t"+"humangenome_sites"+"\t"+"virus_genome"+"\t"+"virus_sites"+"\t"+"support_reads_num"+"\t"+"support_reads_contain"+"\n")
        f.write(sites_cluster_result)




if __name__ == "__main__":
    file = sys.argv[1]
    results =sites_info(file)
    #print(results)
    cluster_human = sites_cluster("blast_sortpy_result")
    last_sum = total_sum('cluster_sort_bycell_setindex')

