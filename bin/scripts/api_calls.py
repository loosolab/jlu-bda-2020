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
    """
    This function creates a dictionary containing the structure of the tree object that angular needs to display it. 
    All analyzed biosources and tfs are in this object. The return of this function is a JSON object.
    """
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'result.csv'))
    results = repository.Repository().read_csv(filename)
    data_ls = []
    data_dict = {}
    
    #create dict containing the biosources and tfs
    for index in results.index:
        biosource = results["biosource"][index] +"_"+ results["genome"][index]
        tf = results["tf"][index] + "_" + results["chr"][index]
        if biosource not in data_dict:
            data_dict[biosource]=[]
        if tf not in data_dict[biosource]:
            data_dict[biosource].append(tf)
    
    #create tree object with the data of the dict
    for biosource in data_dict:
        tree_node = {"item": biosource, "type":"", "belongsTo":"", "checked": False, "children":[]}
        for tf in data_dict[biosource]:
            inner_tree_node = {"item": tf, "type":"tf", "belongsTo":biosource, "checked": False, "children":[]}
            tree_node["children"].append(inner_tree_node)
        data_ls.append(tree_node)
    
    return {"data": data_ls}


def getChecked(data):
    """
    This function gets the modified tree object back from the visualization and analyzes what biosources and tfs are selected by the user.
    The return of this function is a dictionary containing the selected biosources and tfs.
    """
    whats_checked_bio_tf={}
    #analyze what was checked by the user, all biosources only checked tfs
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
    
    #remove empty biosources
    only_checked={}
    for biosource in whats_checked_bio_tf:
        if whats_checked_bio_tf[biosource]!=[]:
            only_checked[biosource]=whats_checked_bio_tf[biosource]
            
    return only_checked
    

def getRawData(checked_data):
    """
    This function creates a dictionary containing the biosources, tfs, both CHIP-seq and ATAC-seq scores and the weights from the analysis part for these.
    This dictionary gets returned as a JSON object.
    """
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'results', 'result.csv'))
    results = repository.Repository().read_csv(filename)
    rawdata_dict = {}
    
    for biosource in checked_data:
        rawdata_dict[biosource]={}
        
        # extract path, weights and means for each tf from the result.csv
        for tf in checked_data[biosource]:
            path_for_tf = results.loc[(results["tf"] == tf.split(sep="_")[0]) & (results["chr"] == tf.split(sep="_")[1]) & (results["biosource"] == biosource.split(sep="_")[0]) & (results["genome"] == biosource.split(sep="_")[1])]["path"].iloc[0]
            weights = results.loc[(results["tf"] == tf.split(sep="_")[0]) & (results["chr"] == tf.split(sep="_")[1]) & (results["biosource"] == biosource.split(sep="_")[0]) & (results["genome"] == biosource.split(sep="_")[1])]["weights"].tolist()
            means = results.loc[(results["tf"] == tf.split(sep="_")[0]) & (results["chr"] == tf.split(sep="_")[1]) & (results["biosource"] == biosource.split(sep="_")[0]) & (results["genome"] == biosource.split(sep="_")[1])]["means"].tolist()
            new_weights=[]
            
            # add all weight into a list and round them
            for weight in weights:
                new_weights.append(round(weight,5))
            
            atac_ls= []
            chip_ls=[]
            
            # add means into a list and round them
            for mean in means:
                atac_ls.append(round(mean[0],5))
                chip_ls.append(round(mean[1],5))
                 
            #secure that path exists
            if os.path.exists(path_for_tf):
                csvfile = open(os.path.join(path_for_tf,tf.split(sep="_")[0]+".csv")) 
                data = list(csv.reader(csvfile, delimiter=","))
                csvfile.close()
                x = []
                y = []
                
                # add CHIP-seq and ATAC-seq scores as x and y coordinates to lists
                for row in data:
                    x.append(round(float(row[0]),3))
                    y.append(round(float(row[1]),3))
                
                # add all data for each tf into the dictionary
                rawdata_dict[biosource][tf]={}
                rawdata_dict[biosource][tf]["rawData"]=[]
                rawdata_dict[biosource][tf]["analysisData"]=[]
                rawdata_dict[biosource][tf]["analysisData"].append(atac_ls)
                rawdata_dict[biosource][tf]["analysisData"].append(chip_ls)
                rawdata_dict[biosource][tf]["analysisData"].append(new_weights)
                rawdata_dict[biosource][tf]["rawData"].append([x, y])


    return {"data": rawdata_dict}
