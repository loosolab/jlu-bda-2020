# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 16:53:09 2021

@author: Spyro
"""

import repository
import pandas as pd
import os
import csv
# =============================================================================
# def parse_result_csv():
#     filename = "result.csv"
#     
#     biosource_ls=[]
#     data_ls=[]
#     all_tf_ls=[]
#     results = repository.repository().read_csv(filename)
#    # print(results)
#    
#     
#     for index in results.index:
#         if results["biosource"][index] not in biosource_ls:
#             if biosource_ls != []:
#                 biosource_ls.append(all_tf_ls)
#                 data_ls.append(biosource_ls)
#             biosource_ls=[]
#             biosource_ls.append(results["biosource"][index])
#         tf_path_ls = []
#         tf_path_ls.append(results["tf"][index])
#         tf_path_ls.append(str(index+1))
#         all_tf_ls.append(tf_path_ls)
#         #print(results["biosource"][index], results["tf"][index])
#     biosource_ls.append(all_tf_ls)
#     data_ls.append(biosource_ls)
#     full_dict = {"data": data_ls}
#     return full_dict
# =============================================================================

# =============================================================================
# def get_biosource_list_for_tree():
#     filename = "result.csv"
#     results = repository.repository().read_csv(filename)
#     data_ls = []
#     tree_node = {"item": "", "type":"", "belongsTo":"", "checked": False, "children":[]}
#     inner_tree_node = {"item": "", "type":"", "belongsTo":"", "checked": False, "children":[]}
#     for index in results.index:
#         if results["biosource"][index] not in tree_node["item"]:
#             if tree_node["item"]!= "":
#                 data_ls.append(tree_node)
#             tree_node["item"] = results["biosource"][index]
#         inner_tree_node["item"] = results["tf"][index]
#         inner_tree_node["type"] = "tf"
#         inner_tree_node["belongsTo"] = results["biosource"][index]
#         print(inner_tree_node)
#         tree_node["children"].append(inner_tree_node)
#         
#         
#     data_ls.append(tree_node)
#     #print()
#     #print(tree_node)
#     #print(data_ls)
#     
#     return {"data": data_ls}
# =============================================================================



def get_biosource_list_for_tree():
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'result.csv'))
    results = repository.Repository().read_csv(filename)
    data_ls = []
    data_dict = {}
    for index in results.index:
        biosource = results["biosource"][index]
        tf = results["tf"][index]
        if biosource not in data_dict:
            data_dict[biosource]=[]
        if tf not in data_dict[biosource]:
            data_dict[biosource].append(tf)
    #print(data_dict)
    for biosource in data_dict:
        tree_node = {"item": biosource, "type":"", "belongsTo":"", "checked": False, "children":[]}
        for tf in data_dict[biosource]:
            inner_tree_node = {"item": tf, "type":"tf", "belongsTo":biosource, "checked": False, "children":[]}
            tree_node["children"].append(inner_tree_node)
        data_ls.append(tree_node)
    return {"data": data_ls}

#mydata =  [{'item': 'GM12878', 'type': '', 'belongsTo': '', 'checked': False, 'children': [{'item': 'ARID3A_ENCFF003VDB', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': True, 'children': []}, {'item': 'ARNT_ENCFF758RQJ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': True, 'children': []}]}]

def getChecked(data):

    #print(data)
    whats_checked_bio_tf={}
    for biosource_obj in data:
        biosource = biosource_obj["item"]
        whats_checked_bio_tf[biosource]=[]
        if biosource_obj["checked"]:
            for tf_obj in biosource_obj["children"]:
                tf = tf_obj["item"]
                whats_checked_bio_tf[biosource].append(tf)
        else:
            for tf_obj in biosource_obj["children"]:
                if tf_obj["checked"]:
                    whats_checked_bio_tf[biosource].append(tf_obj["item"])
    #print(whats_checked_bio_tf)
    return whats_checked_bio_tf
    

# =============================================================================
#     for biosource_obj in data:
#         #print(biosource_obj)
#         biosource = biosource_obj["item"]
#         #print(biosource)
# 
#         data_dict["biosource"]= biosource
#         data_dict["tf_list"]=[]
#         for tf_obj in biosource_obj["children"]:
#             if tf_obj["checked"]:
#                 tempdict = {}
#                 tf = tf_obj["item"]
#                 tempdict["tf"]= tf
#                 tempdict["path"] = "path aus csv"
#                 data_dict["tf_list"].append(tempdict)
#         data_ls.append(data_dict)
# =============================================================================
    #print(data_dict)
    # IMPORTANT!! Stay or not stay?
# =============================================================================
# def getGraphs(data):
#     path = "F:\\Uni\\visualization-app\\src\\assets\\img\\plots"
#     path_dict = {}
#     for biosource in data:
#         path_dict[biosource]={}
#         #biosource if missing folder structure
#         for tf in data[biosource]:
#             path_dict[biosource][tf]=os.listdir(path+"\\"+tf)
# 
#     return {"data": path_dict}
# =============================================================================

def getRawData(checked_data):
    #print(checked_data)
    #filepath
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'result.csv'))
    results = repository.Repository().read_csv(filename)
    rawdata_dict = {}
    
    my_local_path = "results\\Genome\\plots\\"
    #print(results.keys())
    
    for biosource in checked_data:
        rawdata_dict[biosource]={}
        for tf in checked_data[biosource]:
            path_for_tf = results.loc[(results["tf"] == tf) & (results["biosource"] == biosource)]["path"].iloc[0]
            weights = results.loc[(results["tf"] == tf) & (results["biosource"] == biosource)]["weights"].tolist()
            means = results.loc[(results["tf"] == tf) & (results["biosource"] == biosource)]["means"].tolist()
            #covariances = results.loc[(results["tf"] == tf) & (results["biosource"] == biosource)]["covariances"]
            print(weights)
            #print(means)
            new_weights=[]
            for weight in weights:
                new_weights.append(round(weight,5))
            analysisData_ls = []
            #print("start")
            atac_ls= []
            chip_ls=[]
            for mean in means:
                #print( means[index])
                
                atac_ls.append(round(mean[0],5))
                chip_ls.append(round(mean[1],5))
                # column_ls --> weight, atac, chip
                #analysisData_ls.append([weights[index], means[index][0], means[index][1]])
                #analysisData_ls.append(column_ls)
            print(atac_ls)
            print(chip_ls)
            
            
            print(analysisData_ls)
            #print(means)
            #print(covariances)
            
            #### needs to be changed from local to full           
            #path_for_tf= path_for_tf.split("/")[-1]
            print(path_for_tf)
            #secure that path exists
            if os.path.exists(path_for_tf):
                csvfile = open(os.path.join(path_for_tf,tf+".csv")) 
                data = list(csv.reader(csvfile, delimiter=","))
                csvfile.close()
                x = []
                y = []
                #z = []
                for row in data:
                    ### Setting Max Value to 100 --> only for testing!!
                    #if float(row[0]) <= 100 and float(row[1])<=100:
                    x.append(round(float(row[0]),3))
                    y.append(round(float(row[1]),3))
                        #z.append(round(float(row[2]),3))
                rawdata_dict[biosource][tf]={}
                rawdata_dict[biosource][tf]["rawData"]=[]
                rawdata_dict[biosource][tf]["analysisData"]=[]
                rawdata_dict[biosource][tf]["analysisData"].append(atac_ls)
                rawdata_dict[biosource][tf]["analysisData"].append(chip_ls)
                rawdata_dict[biosource][tf]["analysisData"].append(new_weights)
                
                #rawdata_dict[biosource][tf]["analysisData"].append([weights, means, covariances])
                rawdata_dict[biosource][tf]["rawData"].append([x, y])
                #print(len(x), len(y), len(z))
            
    #print(rawdata_dict )
    return {"data": rawdata_dict}

#getRawData(getChecked(mydata))