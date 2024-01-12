## Title     : The way for extract integration
## Objective : Result of Human-chr Human-site Virus-type Virus-breakpoint
## Created by: fpeng
## Created on: 2024/01/12


import sys

def only_bk(path):
    with open(path,"r") as f:
        f1 = f.readlines()
    result = ""
    num = 1
    last_content = f1[0]
    for i in f1[1:]:

        if i.split("\t")[0] == last_content.split("\t")[0] and abs(int(i.split("\t")[1])-int(last_content.split("\t")[1]))<=20:
            num +=1
        elif i==f1[-1] and i.split("\t")[0] == last_content.split("\t")[0] and abs(int(i.split("\t")[1])-int(last_content.split("\t")[1]))<=20:
            result += last_content.split("\t")[0]+"\t" + last_content.split("\t")[1]+"\t"+ str(num+1)+"\n"
        elif i==f1[-1] and i.split("\t")[0] != last_content.split("\t")[0]:
            result += last_content.split("\t")[0]+"\t" + last_content.split("\t")[1]+"\t"+ str(num)+"\n"+i.split("\t")[0] + "\t" + i.split("\t")[1] + "\t" + str(num) + "\n"
        else:
            result += last_content.split("\t")[0]+"\t" + last_content.split("\t")[1]+"\t"+ str(num)+"\n"
            num =1
            last_content = i
    print(result)


if __name__ == "__main__":
    argv = sys.argv[1]
    last_result = only_bk(argv)
