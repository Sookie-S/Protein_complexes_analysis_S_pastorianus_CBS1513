#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 20:11:42 2020

@author: Soukaina Timouma
"""


#%%
import os
import re
import matplotlib.pyplot as plt
import numpy as np
#%%

os.chdir("/home/sookie/Documents/RNAseq_analysis/9_Young_paralogs")

#%%
os.getcwd()



#%%
#
#
#
#    with open("results/homodimers.csv","w") as out:
#    out.write("Number	Complex accession	Recommended name	ID	S. eubayanus-like allele	S.cerevisiae-like allele	S. eubayanus absolute expression at 22 degrees	S. cerevisiae absolute expression at 22 degrees	Conclusion	Comments\n")
#                
#                out.write("\t".join(dico_comp_expr[comp][0])+"\t"+str(conc)+"\t"+str(comment)+"\n")
#
#
#
#with open("results/heterodimers.csv","w") as out:
#    out.write("Number	Complex accession	Recommended name	ID	S. eubayanus-like allele	S.cerevisiae-like allele	S. eubayanus absolute expression at 22 degrees	S. cerevisiae absolute expression at 22 degrees	Conclusion	Comments	Hypothesis\n")
#
#                out.write("\t".join(dico_comp_expr[comp][0])+"\t"+str(conc)+"\t"+str(comment0)+"\t"+str(hyp)+"\n")
#                out.write("\t".join(dico_comp_expr[comp][1])+"\t"+str(conc)+"\t"+str(comment1)+"\t"+str(hyp)+"\n") 
#
#
#    with open("results/trimers.csv","w") as out:
#        out.write("Number	Complex accession	Recommended name	ID	S. eubayanus-like allele	S.cerevisiae-like allele	S. eubayanus absolute expression at 22 degrees	S. cerevisiae absolute expression at 22 degrees	Conclusion	Comments	Hypothesis\n")
#
#
#                out.write("\t".join(dico_comp_expr[comp][0])+"\t"+str(conc)+"\t"+str(comment0)+"\t"+str(hyp)+"\n")
#                out.write("\t".join(dico_comp_expr[comp][1])+"\t"+str(conc)+"\t"+str(comment1)+"\t"+str(hyp)+"\n") 
#                out.write("\t".join(dico_comp_expr[comp][2])+"\t"+str(conc)+"\t"+str(comment2)+"\t"+str(hyp)+"\n") 

#############################################################################################################

#%%absolute expression


dico_comp_expr = {}
with open("results/Table_complexes_in_cbs1513_FINAL.csv","r") as exp:
    next(exp)
    for line in exp:
        line = line.split("\n")[0]
        line = line.split(",")
        
        if line[0] not in dico_comp_expr.keys():
            x=0
            dico_comp_expr[line[0]] = {x:line}
        else:
            x+=1
            dico_comp_expr[line[0]][x] = line
#        print(line)
        
# 
#    
#%%
#print(dico_comp_expr[comp][0])
#se_expr
#dico_comp_expr['CPX-1355']    
#dico_absolute_expr_13w
#
#%%              
# 13w, 22w, 30w, 13aa, 22aa, 30aa, 13leu, 22leu, 30leu, 13eth, 22eth, 30eth 
temperatures_columns = []
temperatures_columns.append([1,2,3,"13C in wort","13w"])
temperatures_columns.append([4,5,6,"22C in wort","22w"])
temperatures_columns.append([7,8,9,"30C in wort","30w"])
temperatures_columns.append([10,11,12,"13C in SD media","13aa"])
temperatures_columns.append([13,14,15,"22C in SD media","22aa"])
temperatures_columns.append([16,17,18,"30C in SD media","30aa"])
temperatures_columns.append([19,20,21,"13C in SD media w/o leucine","13leu"])
temperatures_columns.append([22,23,24,"22C in SD media w/o leucine","22leu"])
temperatures_columns.append([25,26,27,"30C in SD media w/o leucine","30leu"])
temperatures_columns.append([28,29,30,"13C in SD media + 6% ethanol","13eth"])
temperatures_columns.append([31,32,33,"22C in SD media + 6% ethanol","22eth"])
temperatures_columns.append([34,35,36,"30C in SD media + 6% ethanol","30eth"])

threshold = 3
print("######################### THRESHOLD :",threshold)

#
##%%
#print(comp)
#print("Seub: ",dico_comp_expr[comp][0][4])
#print(globals()['dico_absolute_expr_'+str(temperature[4])][dico_comp_expr[comp][0][4]])
##print("Scer: ",dico_comp_expr[comp][0][5])
##print(globals()['dico_absolute_expr_'+str(temperature[4])][dico_comp_expr[comp][0][5]])
#
##%%
res_homodimers = {}
res_heterodimers = {}
for temperature in temperatures_columns:
    print("\n\n")
    print("####################### ",temperature[3]," #####################################\n")
#    print(temperature)
    globals()['dico_absolute_expr_'+str(temperature[4])] = {}
    with open("files/annotated_table_expression_all_conditions.csv","r") as annot:
        next(annot)
        for line in annot:
            line = line.split("\n")[0]
            line = line.split(",")
            
            mean_expr = (float(line[temperature[0]])+float(line[temperature[1]])+float(line[temperature[2]]))/3
#            print(line[0],float(line[temperature[0]]),float(line[temperature[1]]),float(line[temperature[2]]))
            sp = line[0].replace("_","")
#            print(mean_expr)
            globals()['dico_absolute_expr_'+str(temperature[4])][sp] = mean_expr          
    annot.close()
    
    frTOfr = 0
    frTOpr = 0
    frTOchi = 0
    frTOuni = 0
    prTOpr = 0
    prTOchi = 0
    prTOuni = 0
    cTOc = 0
    uniTOuni =0
    
    pr = 0
    fr = 0
    chi = 0
    uni = 0
    
    homodimers_issue = []
    count_except_homo = 0
    print("CASE HOMODIMER\n")
    
    
    for comp in dico_comp_expr.keys():
        if len(dico_comp_expr[comp].values()) == 1:
#            print(dico_comp_expr[comp])
            
            se_allele = dico_comp_expr[comp][0][3]
            sc_allele = dico_comp_expr[comp][0][4]
            
            try:
                se_expr = globals()['dico_absolute_expr_'+str(temperature[4])][se_allele]
            except:
                homodimers_issue.append(se_allele)
                count_except_homo+=1
                se_expr = 0
                pass
            try:
                sc_expr = globals()['dico_absolute_expr_'+str(temperature[4])][sc_allele]
            except:
                homodimers_issue.append(sc_allele)
                count_except_homo+=1
                sc_expr = 0
                pass 
            
            if se_allele != "NA" and sc_allele != "NA" :
                conc = "Fully redundant"
                fr+=1
                
                if comp not in res_homodimers.keys():
                    res_homodimers[comp] = {'genomic': conc}
                
                if se_expr >= threshold*float(sc_expr):
                    comment = "S. eubayanus +++"
                    frTOuni+=1
                    res_homodimers[comp][str(temperature[4])] = 'Unispecific Se'
                elif float(sc_expr) >= threshold*float(se_expr):
                    comment = "S. cerevisiae +++"
                    res_homodimers[comp][str(temperature[4])] = 'Unispecific Sc'
                    frTOuni+=1
                else:
                    res_homodimers[comp][str(temperature[4])] = 'Fully redundant'
                    comment = "Both alleles"
                    frTOfr+=1
            else:
                conc = "Unispecific"
                if comp not in res_homodimers.keys():
                    res_homodimers[comp] = {'genomic': conc}
                if se_allele != "NA":
                    comment = "S. eubayanus"
                    res_homodimers[comp][str(temperature[4])] = 'Unispecific Se'
                    uni+=1
                else:
                    comment = "S. cerevisiae"
                    res_homodimers[comp][str(temperature[4])] = 'Unispecific Sc'
                    uni+=1


#                print(str(conc)+"\texpression --> "+str(comment)+"\n")

    tot = fr + uni
    print("FR: ",fr)
    print("FR to FR: ",frTOfr)
    print("FR to UNI: ",frTOuni)
    print("UNI: ",uni)
    print("Total: ",tot)    
    
    print("\n")  
    
       
              
    print("CASE HETERODIMER\n")
    frTOfr = 0
    frTOpr = 0
    frTOchi = 0
    frTOuni = 0
    prTOpr = 0
    prTOchi = 0
    prTOuni = 0
    cTOc = 0
    uniTOuni =0
    
    pr = 0
    fr = 0
    chi = 0
    uni = 0    
    
    
    heterodimers_issue = []
    count_except_hetero = 0
    for comp in dico_comp_expr.keys():

        if len(dico_comp_expr[comp].values()) == 2:
#                print(dico_comp_expr[comp][0])
#                print(dico_comp_expr[comp][1])
            
            
            se_allele1 = dico_comp_expr[comp][0][3]
            sc_allele1 = dico_comp_expr[comp][0][4]
            se_allele2 = dico_comp_expr[comp][1][3]
            sc_allele2 = dico_comp_expr[comp][1][4]   
            
            try:
                se_expr1 = globals()['dico_absolute_expr_'+str(temperature[4])][se_allele1]
            except:
                heterodimers_issue.append(se_allele1)
                count_except_hetero+=1
                se_expr1 = 0
                pass                
            
            try:
                sc_expr1 = globals()['dico_absolute_expr_'+str(temperature[4])][sc_allele1]  
            except:
                heterodimers_issue.append(sc_allele1)
                count_except_hetero+=1
                sc_expr1 = 0
                pass   
            
            try:
                count_except_hetero+=1
                se_expr2 = globals()['dico_absolute_expr_'+str(temperature[4])][se_allele2]
            except:
                heterodimers_issue.append(se_allele2)
                count_except_hetero+=1
                se_expr2 = 0
                pass   

            try:
                sc_expr2 = globals()['dico_absolute_expr_'+str(temperature[4])][sc_allele2]    
            except:
                heterodimers_issue.append(sc_allele2)
                count_except_hetero+=1
                sc_expr2 = 0
                pass   

           
            if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" :
                conc = "Fully redundant"
                fr+=1
                
                if comp not in res_heterodimers.keys():
                    res_heterodimers[comp]= {'genomic': conc} 
                if float(se_expr1) >= threshold*float(sc_expr1):
                    comment0 = "S. eubayanus +++"
                elif float(sc_expr1) >= threshold*float(se_expr1):
                    comment0 = "S. cerevisiae +++"
                else:
                    comment0 = "Both alleles"
                    
                if float(se_expr2) >= threshold*float(sc_expr2):
                    comment1 = "S. eubayanus +++"
                elif float(sc_expr2) >= threshold*float(se_expr2):
                    comment1 = "S. cerevisiae +++"
                else:
                    comment1 = "Both alleles"
                    
                if comment0 == "S. cerevisiae +++" and comment1 == "S. cerevisiae +++":
                    hyp = "Unispecific S. cerevisiae"
                    res_heterodimers[comp][str(temperature[4])] ='Unispecific Sc'
                elif comment0 == "S. eubayanus +++" and comment1 == "S. eubayanus +++":
                    hyp = "Unispecific S. eubayanus"
                    res_heterodimers[comp][str(temperature[4])] ='Unispecific Se'
                    
                elif (comment0 == "S. cerevisiae +++" and comment1 == "S. eubayanus +++"):
                    hyp =  "Chimeric"
                    res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Sc P2-Se"
                elif (comment1 == "S. cerevisiae +++" and comment0 == "S. eubayanus +++") :
                    hyp =  "Chimeric"
                    res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Se P2-Sc"                                       
                    
                elif (comment0 == "S. cerevisiae +++" and comment1 == "Both alleles"):
                    hyp = "Partially redundant - one subunit S. cerevisiae"
                    res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P1-Sc"
                elif (comment0 == "Both alleles" and comment1 == "S. cerevisiae +++"):
                    hyp = "Partially redundant - one subunit S. cerevisiae"
                    res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P2-Sc"
                elif (comment0 == "S. eubayanus +++" and comment1 == "Both alleles"):
                    hyp = "Partially redundant - one subunit S. eubayanus"
                    res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P1-Se"
                elif (comment0 == "Both alleles" and comment1 == "S. eubayanus +++"):
                    hyp = "Partially redundant - one subunit S. eubayanus"
                    res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P2-Se"                    
                    
                elif comment0 == "Both alleles" and comment1 == "Both alleles":
                    hyp = "Fully redundant"
                    res_heterodimers[comp][str(temperature[4])] =hyp
                else:
#                        print(comment0,comment1)
                    print("ISSUE CASE 1")
                    break
                   
            elif (se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 != "NA") or (se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA") or (se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA") or (se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA"):
                conc = "Partially redundant"
                pr+=1 
                
                if comp not in res_heterodimers.keys():
                    res_heterodimers[comp]= {'genomic': conc}                 
                try:
                    if float(se_expr1) >= threshold*float(sc_expr1):
                        comment0 = "S. eubayanus +++"
                    elif float(sc_expr1) >= threshold*float(se_expr1):
                        comment0 = "S. cerevisiae +++"
                    else:
                        comment0 = "Both alleles"
                except:
                    pass
                try:
                    if float(se_expr2) >= threshold*float(sc_expr2):
                        comment1 = "S. eubayanus +++"
                    elif float(sc_expr2) >= threshold*float(se_expr2):
                        comment1 = "S. cerevisiae +++"
                    else:
                        comment1 = "Both alleles"                    
                except:
                    pass
                    
                if se_allele1 != "NA" and sc_allele1 == "NA":
                    comment0 = "S. eubayanus"
                    if comment1 == "S. eubayanus +++":
                        hyp = "Unispecific S. eubayanus"
                        res_heterodimers[comp][str(temperature[4])] ='Unispecific Se'    
                    elif comment1 == "S. cerevisiae +++":
                        hyp = "Chimeric"
                        res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Se P2-Sc"
                    elif comment1 == "Both alleles":
                        hyp = "Partially redundant - one subunit S. eubayanus"
                        res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P1-Se"
                    else:
#                            print(comment0,comment1)
                        print("ISSUE CASE 2")
                        break
                    
                if se_allele1 == "NA" and sc_allele1 != "NA":
                    comment0 = "S. cerevisiae"
                    if comment1 == "S. cerevisiae +++":
                        hyp = "Unispecific S. cerevisiae"
                        res_heterodimers[comp][str(temperature[4])] ='Unispecific Sc'
                    elif comment1 == "S. eubayanus +++":
                        hyp = "Chimeric"
                        res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Sc P2-Se"
                    elif comment1 == "Both alleles":
                        hyp = "Partially redundant - one subunit S. cerevisiae"
                        res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P1-Sc"
                    else:
#                            print(comment0,comment1)
                        print("ISSUE CASE 2 bis")
                        break
                    
                if se_allele2 != "NA" and sc_allele2 == "NA":
                    comment1 = "S. eubayanus"
                    if comment0 == "S. eubayanus +++":
                        hyp = "Unispecific S. eubayanus"
                        res_heterodimers[comp][str(temperature[4])] ='Unispecific Se'
                    elif comment0 == "S. cerevisiae +++":
                        hyp = "Chimeric"
                        res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Sc P2-Se"
                    elif comment0 == "Both alleles":
                        hyp = "Partially redundant - one subunit S. eubayanus"
                        res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P2-Se"
                    else:
#                            print(comment0,comment0)
                        print("ISSUE CASE 3")
                        break                        
                    
                if se_allele2 == "NA" and sc_allele2 != "NA":
                    comment1 = "S. cerevisiae"
                    if comment0 == "S. cerevisiae +++":
                        hyp = "Unispecific S. cerevisiae"
                        res_heterodimers[comp][str(temperature[4])] ='Unispecific Sc'
                    elif comment0 == "S. eubayanus +++":
                        hyp = "Chimeric"
                        res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Se P2-Sc"
                    elif comment0 == "Both alleles":
                        hyp = "Partially redundant - one subunit S. cerevisiae"
                        res_heterodimers[comp][str(temperature[4])] ="Partially redundant: P2-Sc"
                    else:
#                            print(comment0,comment0)
                        print("ISSUE CASE 4")
                        break
                     
            elif (se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 == "NA" and sc_allele2 != "NA") or (se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA"):
                conc = "Chimeric"   
                hyp = "Chimeric"
                chi+=1
                
                if comp not in res_heterodimers.keys():
                    res_heterodimers[comp]= {'genomic': conc} 
                if se_allele1 != "NA":
                    comment0 = "S. eubayanus"
                    res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Se P2-Sc"
                else:
                    comment0 = "S. cerevisiae"
                    res_heterodimers[comp][str(temperature[4])] ="Chimeric: P1-Sc P2-Se"
                if se_allele2 != "NA":
                    comment1 = "S. eubayanus"
                else:
                    comment1 = "S. cerevisiae"
                    
            else:
                uni+=1
                conc = "Unispecific"
                
                if comp not in res_heterodimers.keys():
                    res_heterodimers[comp]= {'genomic': conc} 
                if se_allele1 != "NA":
                    comment0 = "S. eubayanus"
                    hyp = "Unispecific S. eubayanus"
                    res_heterodimers[comp][str(temperature[4])] ='Unispecific Se'
                else:
                    comment0 = "S. cerevisiae"
                    hyp = "Unispecific S. cerevisiae"
                    res_heterodimers[comp][str(temperature[4])] ='Unispecific Sc'
                if se_allele2 != "NA":
                    comment1 = "S. eubayanus"
                    hyp = "Unispecific S. eubayanus"
                else:
                    comment1 = "S. cerevisiae"
                    hyp = "Unispecific S. cerevisiae"
               
            
            if conc == "Fully redundant" and hyp == "Fully redundant":
                frTOfr +=1      
            elif conc == "Fully redundant" and "Partially redundant" in hyp:
                frTOpr +=1                       
            elif conc == "Fully redundant" and hyp == "Chimeric":
                frTOchi +=1
            elif conc == "Fully redundant" and "Unispecific" in hyp:
                frTOuni +=1
            elif conc == "Partially redundant" and "Partially redundant" in hyp:
                prTOpr +=1                
            elif conc == "Partially redundant" and hyp == "Chimeric":
                prTOchi+=1
            elif conc == "Partially redundant" and "Unispecific" in hyp:
                prTOuni+=1
            elif conc == "Chimeric" and hyp == "Chimeric":
                cTOc+=1
            elif conc == "Unispecific"  and "Unispecific" in hyp:
                uniTOuni+=1
            else:
                print("OTHER CASE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
                break
            
    tot = fr + pr + uni+ chi
    print("FR: ",fr)
    print("FR to FR: ",frTOfr)
    print("FR to PR: ",frTOpr)
    print("FR to UNI: ",frTOuni)
    print("FR to CHI: ",frTOchi)
    print("PR: ",pr)
    print("PR to PR: ",prTOpr)
    print("PR to CHI: ",prTOchi)
    print("PR to UNI: ",prTOuni)
    print("CHI: ",chi)
    print("CHI to CHI: ",cTOc)
    print("UNI: ",uni)
    print("UNI to UNI: ",uniTOuni)
    print("Total: ",tot)
    
    if fr != frTOfr+frTOpr+frTOuni+frTOchi:
        print("ERROR COUNT FR")
        break
    
    if pr != prTOpr+prTOchi+prTOuni:
        print("ERROR COUNT PR")
        break        

    
    


    print("\n\nCASE TRIMERS\n")        
    ### TRIMERS
    frTOfr = 0
    frTOpr = 0
    frTOchi = 0
    frTOuni = 0
    prTOpr = 0
    prTOchi = 0
    prTOuni = 0
    cTOc = 0
    uniTOuni =0
    
    pr = 0
    fr = 0
    chi = 0
    uni = 0

    trimers_issue = []
    count_except_tri = 0
    for comp in dico_comp_expr.keys():

        if len(dico_comp_expr[comp].values()) == 3:
#                print(dico_comp_expr[comp][0])
#                print(dico_comp_expr[comp][1])
#                print(dico_comp_expr[comp][2])
            
            
            se_allele1 = dico_comp_expr[comp][0][3]
            sc_allele1 = dico_comp_expr[comp][0][4]
            se_allele2 = dico_comp_expr[comp][1][3]
            sc_allele2 = dico_comp_expr[comp][1][4]       
            se_allele3 = dico_comp_expr[comp][2][3]
            sc_allele3 = dico_comp_expr[comp][2][4] 
            
            try:
                se_expr1 = globals()['dico_absolute_expr_'+str(temperature[4])][se_allele1]
            except:
                trimers_issue.append(se_allele1)
                count_except_tri+=1
                se_expr1 = 0
                pass
            try:
                sc_expr1 = globals()['dico_absolute_expr_'+str(temperature[4])][sc_allele1]   
            except:
                trimers_issue.append(sc_allele1)
                count_except_tri+=1
                sc_expr1 = 0
                pass
            try:
                se_expr2 = globals()['dico_absolute_expr_'+str(temperature[4])][se_allele2]
            except:
                trimers_issue.append(se_allele2)
                count_except_tri+=1
                se_expr2 = 0
                pass
            try:
                sc_expr2 = globals()['dico_absolute_expr_'+str(temperature[4])][sc_allele2]  
            except:
                trimers_issue.append(sc_allele2)
                count_except_tri+=1
                sc_expr2 = 0
                pass
            try:
                se_expr3 = globals()['dico_absolute_expr_'+str(temperature[4])][se_allele3]
            except:
                trimers_issue.append(se_allele3)
                count_except_tri+=1
                se_expr3 = 0
                pass
            try:
                sc_expr3 = globals()['dico_absolute_expr_'+str(temperature[4])][sc_allele3]     
            except:
                trimers_issue.append(sc_allele3)
                count_except_tri+=1
                sc_expr3 = 0
                pass
            
           
            ############################ FULLY REDUNDANT
            if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 != "NA"  :
                conc = "Fully redundant"
                fr+=1
                if float(se_expr1) >= threshold*float(sc_expr1):
                    comment0 = "S. eubayanus +++"
                elif float(sc_expr1) >= threshold*float(se_expr1):
                    comment0 = "S. cerevisiae +++"
                else:
                    comment0 = "Both alleles"
                    
                if float(se_expr2) >= threshold*float(sc_expr2):
                    comment1 = "S. eubayanus +++"
                elif float(sc_expr2) >= threshold*float(se_expr2):
                    comment1 = "S. cerevisiae +++"
                else:
                    comment1 = "Both alleles"
                    
                if float(se_expr3) >= threshold*float(sc_expr3):
                    comment2 = "S. eubayanus +++"
                elif float(sc_expr3) >= threshold*float(se_expr3):
                    comment2 = "S. cerevisiae +++"
                else:
                    comment2 = "Both alleles"                        

                    
                    
                if comment0 == "S. cerevisiae +++" and comment1 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++" :
                    hyp = "Unispecific S. cerevisiae"
                elif comment0 == "S. eubayanus +++" and comment1 == "S. eubayanus +++" and comment2 == "S. eubayanus +++":
                    hyp = "Unispecific S. eubayanus"
                    
                elif (comment0 == "S. cerevisiae +++" and comment1 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++") or (comment0 == "S. cerevisiae +++" and comment1 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++") or (comment0 == "S. eubayanus +++" and comment1 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++") or (comment0 == "S. cerevisiae +++" and comment1 == "S. eubayanus +++" and comment2 == "S. eubayanus +++") or (comment0 == "S. eubayanus +++" and comment1 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++") or (comment0 == "S. eubayanus +++" and comment1 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++"):
                    hyp =  "Chimeric"
                 
                elif (comment0 == "S. cerevisiae +++" and comment1 == "S. cerevisiae +++" and comment2 == "Both alleles") or (comment0 == "S. cerevisiae +++" and comment1 == "Both alleles" and comment2 == "S. cerevisiae +++") or (comment0 == "Both alleles" and comment1 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++"):
                    hyp = "Partially redundant - two subunits S. cerevisiae"
                elif (comment0 == "S. cerevisiae +++" and comment1 == "Both alleles" and comment2 == "Both alleles") or (comment0 == "Both alleles" and comment1 == "Both alleles" and comment2 == "S. cerevisiae +++") or (comment0 == "Both alleles" and comment1 == "S. cerevisiae +++" and comment2 == "Both alleles"):
                    hyp = "Partially redundant - one subunit S. cerevisiae"                        
                    
                elif (comment0 == "S. eubayanus +++" and comment1 == "S. eubayanus +++" and comment2 == "Both alleles") or (comment0 == "S. eubayanus +++" and comment1 == "Both alleles" and comment2 == "S. eubayanus +++") or (comment0 == "Both alleles" and comment1 == "S. eubayanus +++" and comment2 == "S. eubayanus +++"):
                    hyp = "Partially redundant - two subunits S. eubayanus"
                    
                elif (comment0 == "S. eubayanus +++" and comment1 == "Both alleles" and comment2 == "Both alleles") or (comment0 == "Both alleles" and comment1 == "Both alleles" and comment2 == "S. eubayanus +++") or (comment0 == "Both alleles" and comment1 == "S. eubayanus +++" and comment2 == "Both alleles"):
                    hyp = "Partially redundant - one subunit S. eubayanus"                        
                                            
                elif (comment0 == "S. eubayanus +++" and comment1 == "Both alleles" and comment2 == "S. cerevisiae +++") or (comment0 == "S. eubayanus +++" and comment1 == "S. cerevisiae +++" and comment2 == "Both alleles") or (comment0 == "S. cerevisiae +++" and comment1 == "Both alleles" and comment2 == "S. eubayanus +++") or (comment0 == "S. cerevisiae +++" and comment1 == "S. eubayanus +++" and comment2 == "Both alleles") or (comment0 == "Both alleles" and comment1 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++") or (comment0 == "Both alleles" and comment1 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++"):
                    hyp = "Partially redundant - one subunits S. eubayanus - one ubunit S. cerevisiae - one both "
                    
                elif comment0 == "Both alleles" and comment1 == "Both alleles" and comment2 == "Both alleles":
                    hyp = "Fully redundant"
                else:
#                        print(comment0,comment1,comment2,hyp)
                    print("ISSUE CASE 1")
                    

            ################ CHIMERIC      
            
            ## case 2 Se 1 Sc
            elif (se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 == "NA" and sc_allele3 != "NA") :
                comment0 = "S. eubayanus"
                comment1 = "S. eubayanus"
                comment2 = "S. cerevisiae"
                chi+=1
                conc = "Chimeric"   
                hyp = "Chimeric - 2 subunits S. eubayanus - 1 S. cerevisiae"
            elif (se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 == "NA") :
                comment0 = "S. eubayanus"
                comment1 = "S. cerevisiae"
                comment2 = "S. eubayanus"
                chi+=1
                conc = "Chimeric"   
                hyp = "Chimeric - 2 subunits S. eubayanus - 1 S. cerevisiae"
            elif (se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 != "NA" and sc_allele3 == "NA"):
                comment0 = "S. cerevisiae"
                comment1 = "S. eubayanus"
                comment2 = "S. eubayanus"
                chi+=1
                conc = "Chimeric"   
                hyp = "Chimeric - 2 subunits S. eubayanus - 1 S. cerevisiae"
            
            ## case 1 Se 2 Sc
            elif (se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 == "NA") :
                comment0 = "S. cerevisiae"
                comment1 = "S. cerevisiae"
                comment2 = "S. eubayanus"
                chi+=1
                conc = "Chimeric"   
                hyp = "Chimeric - 1 subunits S. eubayanus - 2 S. cerevisiae"
            elif (se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 == "NA" and sc_allele3 != "NA") :
                comment0 = "S. cerevisiae"
                comment1 = "S. eubayanus"
                comment2 = "S. cerevisiae"
                chi+=1
                conc = "Chimeric"   
                hyp = "Chimeric - 1 subunits S. eubayanus - 2 S. cerevisiae"
            elif (se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 == "NA" and sc_allele3 != "NA"):
                comment0 = "S. eubayanus"
                comment1 = "S. cerevisiae"
                comment2 = "S. cerevisiae"
                chi+=1
                conc = "Chimeric"   
                hyp = "Chimeric - 1 subunits S. eubayanus - 2 S. cerevisiae"
 
                    
                    
             ################ UNISPECIFIC   
             # 3 Se
            elif (se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 != "NA" and sc_allele3 == "NA"):
                comment0 = "S. eubayanus"
                comment1 = "S. eubayanus"
                comment2 = "S. eubayanus"
                uni+=1
                conc = "Unispecific"
                hyp = "Unispecific S. eubayanus"
            ## 3 Sc
            elif (se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 == "NA" and sc_allele3 != "NA"):
                comment0 = "S. cerevisiae"
                comment1 = "S. cerevisiae"
                comment2 = "S. cerevisiae"
                uni+=1
                conc = "Unispecific"
                hyp = "Unispecific S. cerevisiae"


            ############################ PARTIALLY REDUNDANT
            else:
                     
                pr+=1
                conc = "Partially redundant"
                                  
                if se_allele1 != "NA" and sc_allele1 != "NA":
                    if float(se_expr1) >= threshold*float(sc_expr1):
                        comment0 = "S. eubayanus +++"
                    elif float(sc_expr1) >= threshold*float(se_expr1):
                        comment0 = "S. cerevisiae +++"
                    else:
                        comment0 = "Both alleles"
                        
                elif se_allele1 == "NA" and sc_allele1 != "NA":
                    comment0 = "S. cerevisiae"
                elif se_allele1 != "NA" and sc_allele1 == "NA":
                    comment0 = "S. eubayanus"
                else:
                    print("Other case PR")
                    break
                    
                        

                if se_allele2 != "NA" and sc_allele2 != "NA":
                    if float(se_expr2) >= threshold*float(sc_expr2):
                        comment1 = "S. eubayanus +++"
                    elif float(sc_expr2) >= threshold*float(se_expr2):
                        comment1 = "S. cerevisiae +++"
                    else:
                        comment1 = "Both alleles"  
                elif se_allele2 == "NA" and sc_allele2 != "NA":
                    comment1 = "S. cerevisiae"
                elif se_allele2 != "NA" and sc_allele2 == "NA":
                    comment1 = "S. eubayanus"
                else:
                    print("Other case PR2")
                    break                            
                        
                        

                if se_allele3 != "NA" and sc_allele3 != "NA":
                    if float(se_expr3) >= threshold*float(sc_expr3):
                        comment2 = "S. eubayanus +++"
                    elif float(sc_expr3) >= threshold*float(se_expr3):
                        comment2 = "S. cerevisiae +++"
                    else:
                        comment2 = "Both alleles"                    
                elif se_allele3 == "NA" and sc_allele3 != "NA":
                    comment2 = "S. cerevisiae"
                elif se_allele3 != "NA" and sc_allele3 == "NA":
                    comment2 = "S. eubayanus"
                else:
                    print("Other case PR3")
                    break  
      
                                  
                 
                ### CASES when 2 subunits have both allele and 1 is unispecific
                
                # case subunits A: Se - B: ScSe - C: SCSe
                if se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment0 = "S. eubayanus"
                    if comment1 == "S. eubayanus +++" and comment2 == "S. eubayanus +++":
                        hyp = "Unispecific S. eubayanus"
                        
                    elif (comment1 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++") or (comment1 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++") or (comment1 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++"):
                        hyp = "Chimeric"
                        
                    elif (comment1 == "Both alleles" and comment2 == "S. cerevisiae +++"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"
                    elif (comment1 == "S. cerevisiae +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"      
                    elif (comment1 == "Both alleles" and comment2 == "S. eubayanus +++"):
                        hyp = "Partially redundant - two subunits S. eubayanus"
                    elif (comment1 == "S. eubayanus +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - two subunits S. eubayanus"                            
                    elif (comment1 == "Both alleles" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - two both"   
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 2")
                        break
                    
                # case subunits A: Sc - B: ScSe - C: SCSe
                if se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment0 = "S. cerevisiae"
                    if comment1 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++":
                        hyp = "Unispecific S. cerevisiae"
                    elif (comment1 == "S. eubayanus +++" and comment2 == "S. eubayanus +++") or (comment1 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++") or (comment1 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++"):
                        hyp = "Chimeric"
                        
                    elif (comment1 == "Both alleles" and comment2 == "S. eubayanus +++"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"
                    elif (comment1 == "S. eubayanus +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"      
                    elif (comment1 == "Both alleles" and comment2 == "S. cerevisiae +++"):
                        hyp = "Partially redundant - two subunits S. cerevisiae"
                    elif (comment1 == "S. cerevisiae +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - two subunits S. cerevisiae"                            
                    elif (comment1 == "Both alleles" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. cerevisiae - two both"   
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 3")
                        break                    
                
                
                # case subunits A: ScSe - B: ScSe - C: Sc
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 == "NA" and sc_allele3 != "NA" :
#                        comment2 = "S. cerevisiae"
                    if comment0 == "S. cerevisiae +++" and comment1 == "S. cerevisiae +++":
                        hyp = "Unispecific S. cerevisiae"
                    elif (comment0 == "S. eubayanus +++" and comment1 == "S. eubayanus +++") or (comment1 == "S. eubayanus +++" and comment0 == "S. cerevisiae +++") or (comment1 == "S. cerevisiae +++" and comment0 == "S. eubayanus +++"):
                        hyp = "Chimeric"   
                    elif (comment1 == "Both alleles" and comment0 == "S. eubayanus +++"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"
                    elif (comment1 == "S. eubayanus +++" and comment0 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"      
                    elif (comment1 == "Both alleles" and comment0 == "S. cerevisiae +++"):
                        hyp = "Partially redundant - two subunits S. cerevisiae"
                    elif (comment1 == "S. cerevisiae +++" and comment0 == "Both alleles"):
                        hyp = "Partially redundant - two subunits S. cerevisiae"                            
                    elif (comment1 == "Both alleles" and comment0 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. cerevisiae - two both"   
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 4")
                        break                    
                                    
                # case subunits A: ScSe - B: ScSe - C: Se
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 == "NA" :
#                        comment2 = "S. eubayanus"
                    if comment0 == "S. eubayanus +++" and comment1 == "S. eubayanus +++":
                        hyp = "Unispecific S. eubayanus"
                    elif (comment0 == "S. cerevisiae +++" and comment1 == "S. cerevisiae +++") or (comment0 == "S. eubayanus +++" and comment1 == "S. cerevisiae +++") or (comment0 == "S. cerevisiae +++" and comment1 == "S. eubayanus +++"):
                        hyp = "Chimeric"
                        
                    elif (comment1 == "Both alleles" and comment0 == "S. cerevisiae +++"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"
                    elif (comment1 == "S. cerevisiae +++" and comment0 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"      
                    elif (comment1 == "Both alleles" and comment0 == "S. eubayanus +++"):
                        hyp = "Partially redundant - two subunits S. eubayanus"
                    elif (comment1 == "S. eubayanus +++" and comment0 == "Both alleles"):
                        hyp = "Partially redundant - two subunits S. eubayanus"                            
                    elif (comment1 == "Both alleles" and comment0 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - two both"   
                    else:
                        
#                            print(comment0,comment1,comment2,conc,hyp)
                        print("\nISSUE CASE 5\n")
#                            print(dico_comp_expr[comp][0])
#                            print(dico_comp_expr[comp][1])
#                            print(dico_comp_expr[comp][2])
                        break                    

                
                # case subunits A: ScSe - B: Se - C: SCSe
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment1 = "S. eubayanus"
                    if comment0 == "S. eubayanus +++" and comment2 == "S. eubayanus +++":
                        hyp = "Unispecific S. eubayanus"
                    elif (comment0 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++") or (comment0 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++") or (comment0 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++"):
                        hyp = "Chimeric"
                        
                    elif (comment0 == "Both alleles" and comment2 == "S. cerevisiae +++"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"
                    elif (comment0 == "S. cerevisiae +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"      
                    elif (comment0 == "Both alleles" and comment2 == "S. eubayanus +++"):
                        hyp = "Partially redundant - two subunits S. eubayanus"
                    elif (comment0 == "S. eubayanus +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - two subunits S. eubayanus"                            
                    elif (comment0 == "Both alleles" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - two both"   
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 6")
                        break                    
                
                # case subunits A: ScSe - B: Sc - C: SCSe
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment1 = "S. cerevisiae"
                    if comment0 == "S. cerevisiae +++" and comment2 == "S. cerevisiae +++":
                        hyp = "Unispecific S. cerevisiae"
                    elif (comment0 == "S. eubayanus +++" and comment2 == "S. eubayanus +++") or (comment0 == "S. eubayanus +++" and comment2 == "S. cerevisiae +++") or (comment0 == "S. cerevisiae +++" and comment2 == "S. eubayanus +++"):
                        hyp = "Chimeric"
                        
                    elif (comment0 == "Both alleles" and comment2 == "S. eubayanus +++"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"
                    elif (comment0 == "S. eubayanus +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. eubayanus - one S. cerevisiae - one both"      
                    elif (comment0 == "Both alleles" and comment2 == "S. cerevisiae +++"):
                        hyp = "Partially redundant - two subunits S. cerevisiae"
                    elif (comment0 == "S. cerevisiae +++" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - two subunits S. cerevisiae"                            
                    elif (comment0 == "Both alleles" and comment2 == "Both alleles"):
                        hyp = "Partially redundant - one subunit S. cerevisiae - two both"   
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 7")
                        break                     
                
                # case subunits A: ScSe - B: Se - C: Se
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 != "NA" and sc_allele3 == "NA" :
#                        comment1 = "S. eubayanus"
#                        comment2 = "S. eubayanus"
                    
                    if comment0 == "S. eubayanus +++" :
                        hyp = "Unispecific S. eubayanus"
                    
                    elif comment0 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment0 == "Both alleles": 
                        hyp = "Partially redundant -two subunits S. eubayanus - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 8")
                        break 
                    
                # case subunits A: ScSe - B: Sc - C: Sc
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 == "NA" and sc_allele3 != "NA" :
#                        comment1 = "S. cerevisiae"
#                        comment2 = "S. cerevisiae"
                    
                    if comment0 == "S. cerevisiae +++" :
                        hyp = "Unispecific S. cerevisiae"
                    
                    elif comment0 == "S. eubayanus +++" :
                        hyp = "Chimeric"
                    elif comment0 == "Both alleles": 
                        hyp = "Partially redundant -two subunits S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 9")
                        break  

                # case subunits A: Sc - B: ScSe - C: Sc
                if se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 == "NA" and sc_allele3 != "NA" :
#                        comment0 = "S. cerevisiae"
#                        comment2 = "S. cerevisiae"
                    
                    if comment1 == "S. cerevisiae +++" :
                        hyp = "Unispecific S. cerevisiae"
                    
                    elif comment1 == "S. eubayanus +++" :
                        hyp = "Chimeric"
                    elif comment1 == "Both alleles": 
                        hyp = "Partially redundant -two subunits S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 10")
                        break  

                # case subunits A: Se - B: ScSe - C: Se
                if se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 == "NA" :
#                        comment0 = "S. eubayanus"
#                        comment2 = "S. eubayanus"
                    
                    if comment1 == "S. eubayanus +++" :
                        hyp = "Unispecific S. eubayanus"
                    
                    elif comment1 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment1 == "Both alleles": 
                        hyp = "Partially redundant -two subunits S. eubayanus - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 11")
#                            break 
                    

                    
                # case subunits A: Sc - B: Sc - C: ScSe
                if se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment1 = "S. cerevisiae"
#                        comment0 = "S. cerevisiae"
                    
                    if comment2 == "S. cerevisiae +++" :
                        hyp = "Unispecific S. cerevisiae"
                    
                    elif comment2 == "S. eubayanus +++" :
                        hyp = "Chimeric"
                    elif comment2 == "Both alleles": 
                        hyp = "Partially redundant -two subunits S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 12")
                        break  

                # case subunits A: Se - B: Se - C: ScSe
                if se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment1 = "S. eubayanus"
#                        comment0 = "S. eubayanus"
                    
                    if comment2 == "S. eubayanus +++" :
                        hyp = "Unispecific S. eubayanus"
                    
                    elif comment2 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment2 == "Both alleles": 
                        hyp = "Partially redundant -two subunits S. eubayanus - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 13")
                        break 


                # case subunits A: ScSe - B: Sc - C: Se
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 == "NA" :
#                        comment1 = "S. cerevisiae"
#                        comment2 = "S. eubayanus"
                    
                    if comment0 == "S. eubayanus +++" or comment0 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment0 == "Both alleles": 
                        hyp = "Partially redundant -one subunit S. eubayanus - one S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 14")
                        break
                    
                # case subunits A: ScSe - B: Se - C: Sc
                if se_allele1 != "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 == "NA" and sc_allele3 != "NA" :
#                        comment2 = "S. cerevisiae"
#                        comment1 = "S. eubayanus"
                    
                    if comment0 == "S. eubayanus +++" or comment0 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment0 == "Both alleles": 
                        hyp = "Partially redundant -one subunit S. eubayanus - one S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 15")
                        break

                # case subunits A: Sc - B: ScSe - C: Se
                if se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 == "NA" :
#                        comment0 = "S. cerevisiae"
#                        comment2 = "S. eubayanus"
                    
                    if comment1 == "S. eubayanus +++" or comment1 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment1 == "Both alleles": 
                        hyp = "Partially redundant -one subunit S. eubayanus - one S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 16")
                        break                       
                        
                # case subunits A: Se - B: ScSe - C: Sc
                if se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 != "NA" and sc_allele2 != "NA" and se_allele3 == "NA" and sc_allele3 != "NA" :
#                        comment2 = "S. cerevisiae"
#                        comment0 = "S. eubayanus"
                    
                    if comment1 == "S. eubayanus +++" or comment1 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment1 == "Both alleles": 
                        hyp = "Partially redundant -one subunit S. eubayanus - one S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 17")
                        break                        

                # case subunits A: Sc - B: Se - C: ScSe
                if se_allele1 == "NA" and sc_allele1 != "NA" and se_allele2 != "NA" and sc_allele2 == "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment0 = "S. cerevisiae"
#                        comment1 = "S. eubayanus"
                    
                    if comment2 == "S. eubayanus +++" or comment2 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment2 == "Both alleles": 
                        hyp = "Partially redundant -one subunit S. eubayanus - one S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 18")
                        break
                  
                # case subunits A: Se - B: Sc - C: ScSe
                if se_allele1 != "NA" and sc_allele1 == "NA" and se_allele2 == "NA" and sc_allele2 != "NA" and se_allele3 != "NA" and sc_allele3 != "NA" :
#                        comment1 = "S. cerevisiae"
#                        comment0 = "S. eubayanus"
                    
                    if comment2 == "S. eubayanus +++" or comment2 == "S. cerevisiae +++" :
                        hyp = "Chimeric"
                    elif comment2 == "Both alleles": 
                        hyp = "Partially redundant -one subunit S. eubayanus - one S. cerevisiae - one both"
 
                    else:
#                            print(comment0,comment1,comment2)
                        print("ISSUE CASE 19")
                        break

 

                
            if conc == "Fully redundant" and "Fully redundant" in hyp:
                frTOfr +=1      
            elif conc == "Fully redundant" and "Partially redundant" in hyp:
                frTOpr +=1                       
            elif conc == "Fully redundant" and "Chimeric" in hyp:
                frTOchi +=1
            elif conc == "Fully redundant" and "Unispecific" in hyp:
                frTOuni +=1
            elif conc == "Partially redundant" and "Partially redundant" in hyp:
                prTOpr +=1                
            elif conc == "Partially redundant" and "Chimeric" in hyp:
                prTOchi+=1
            elif conc == "Partially redundant" and "Unispecific" in hyp:
                prTOuni+=1
            elif conc == "Chimeric" and "Chimeric" in hyp:
                cTOc+=1
            elif conc == "Unispecific"  and "Unispecific" in hyp:
                uniTOuni+=1
            else:
                print("OTHER CASE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1")
                break
 
    tot = fr + pr + uni+ chi
    print("FR: ",fr)
    print("FR to FR: ",frTOfr)
    print("FR to PR: ",frTOpr)
    print("FR to UNI: ",frTOuni)
    print("FR to CHI: ",frTOchi)
    print("PR: ",pr)
    print("PR to PR: ",prTOpr)
    print("PR to CHI: ",prTOchi)
    print("PR to UNI: ",prTOuni)
    print("CHI: ",chi)
    print("CHI to CHI: ",cTOc)
    print("UNI: ",uni)
    print("UNI to UNI: ",uniTOuni)
    print("Total: ",tot)
    if fr != frTOfr+frTOpr+frTOuni+frTOchi:
        print("ERROR COUNT FR")
        break
    
    if pr != prTOpr+prTOchi+prTOuni:
        print("ERROR COUNT PR")
        break  
    
    
    
#    break
    
        

#%%
res_heterodimers

with open ("results/threshhold_3/Homodimers_behaviour_across_all_conditions.csv","w") as out:
    out.write("Complex\tGenenomic case\t13C in wort\t22C in wort\t30C in wort\t 13C in SD media\t22C in SD media\t30C in SD media\t13C in SD media w/o leucine\t22C in SD media w/o leucine\t30C in SD media w/o leucine\t13C in SD media + 6% ethanol\t22C in SD media + 6% ethanol\t30C in SD media + 6% ethanol\n")
    for comp in res_homodimers:
        out.write(str(comp)+"\t")
        for key, value in res_homodimers[comp].items():
            print(x)
            out.write(str(value)+"\t")
        out.write("\n")
        
        
with open ("results/threshhold_3/Fully_Redundant_heterodimers_behaviour_across_all_conditions.csv","w") as out:
    out.write("Complex\tGenenomic case\t13C in wort\t22C in wort\t30C in wort\t 13C in SD media\t22C in SD media\t30C in SD media\t13C in SD media w/o leucine\t22C in SD media w/o leucine\t30C in SD media w/o leucine\t13C in SD media + 6% ethanol\t22C in SD media + 6% ethanol\t30C in SD media + 6% ethanol\n")
        
    for comp in res_heterodimers:
        if res_heterodimers[comp]['genomic'] == "Fully redundant":
            out.write(str(comp)+"\t")
            for key, value in res_heterodimers[comp].items():
                print(x)
                out.write(str(value)+"\t")
            out.write("\n")
#        break

with open ("results/threshhold_3/Partially_Redundant_heterodimers_behaviour_across_all_conditions.csv","w") as out:
    out.write("Complex\tGenenomic case\t13C in wort\t22C in wort\t30C in wort\t 13C in SD media\t22C in SD media\t30C in SD media\t13C in SD media w/o leucine\t22C in SD media w/o leucine\t30C in SD media w/o leucine\t13C in SD media + 6% ethanol\t22C in SD media + 6% ethanol\t30C in SD media + 6% ethanol\n")
    for comp in res_heterodimers:
        if res_heterodimers[comp]['genomic'] == "Partially redundant":
            out.write(str(comp)+"\t")
            for key, value in res_heterodimers[comp].items():
                print(x)
                out.write(str(value)+"\t")
            out.write("\n")
#        break
        
        

#%%
print(count_except_homo)
print(count_except_hetero)
print(count_except_tri)

###########################################################################################################
#%% 
                   
print(homodimers_issue)
                   
print(heterodimers_issue)

print(trimers_issue)