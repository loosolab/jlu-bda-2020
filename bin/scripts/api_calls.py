# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 16:53:09 2021

@author: Spyro
"""

import repository
import pandas as pd
import os
import csv

def get_biosource_list_for_tree():
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'result.csv'))
    results = repository.Repository().read_csv(filename)
    data_ls = []
    data_dict = {}
    for index in results.index:
        biosource = results["biosource"][index] +"_"+ results["genome"][index]
        tf = results["tf"][index] + "_" + results["chr"][index]
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


def getChecked(data):
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
    return whats_checked_bio_tf
    

def getRawData(checked_data):

    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'result.csv'))
    results = repository.Repository().read_csv(filename)
    rawdata_dict = {}
    #my_local_path = "results\\Genome\\plots\\"
    
    for biosource in checked_data:
        rawdata_dict[biosource]={}
        for tf in checked_data[biosource]:
            path_for_tf = results.loc[(results["tf"] == tf.split(sep="_")[0]) & (results["chr"] == tf.split(sep="_")[1]) & (results["biosource"] == biosource.split(sep="_")[0]) & (results["genome"] == biosource.split(sep="_")[1])]["path"].iloc[0]
            weights = results.loc[(results["tf"] == tf.split(sep="_")[0]) & (results["chr"] == tf.split(sep="_")[1]) & (results["biosource"] == biosource.split(sep="_")[0]) & (results["genome"] == biosource.split(sep="_")[1])]["weights"].tolist()
            means = results.loc[(results["tf"] == tf.split(sep="_")[0]) & (results["chr"] == tf.split(sep="_")[1]) & (results["biosource"] == biosource.split(sep="_")[0]) & (results["genome"] == biosource.split(sep="_")[1])]["means"].tolist()
            #covariances = results.loc[(results["tf"] == tf) & (results["biosource"] == biosource)]["covariances"]
            #print(means)
            new_weights=[]
            for weight in weights:
                new_weights.append(round(weight,5))
            #print("start")
            atac_ls= []
            chip_ls=[]
            for mean in means:
                #print( means[index])
                
                atac_ls.append(round(mean[0],5))
                chip_ls.append(round(mean[1],5))
                # column_ls --> weight, atac, chip
            
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
                    x.append(round(float(row[0]),3))
                    y.append(round(float(row[1]),3))

                rawdata_dict[biosource][tf]={}
                rawdata_dict[biosource][tf]["rawData"]=[]
                rawdata_dict[biosource][tf]["analysisData"]=[]
                rawdata_dict[biosource][tf]["analysisData"].append(atac_ls)
                rawdata_dict[biosource][tf]["analysisData"].append(chip_ls)
                rawdata_dict[biosource][tf]["analysisData"].append(new_weights)
                
                #rawdata_dict[biosource][tf]["analysisData"].append([weights, means, covariances])
                rawdata_dict[biosource][tf]["rawData"].append([x, y])


    return {"data": rawdata_dict}
