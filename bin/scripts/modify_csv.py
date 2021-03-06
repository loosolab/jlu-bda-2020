#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 17:56:52 2021

@author: jan
"""
from scripts.repository import Repository as rp
import pandas as pd
import os

class modifyCSV:
    
    def __init__(self, path):
        """
        

        Parameters
        ----------
        path : path to results directory

        Returns
        -------
        None.

        """
        self.path = path 
    
    def compare(self, parameters):
        """
        This Method modifys an existing result.csv by comparing arguments row by row.
        If a previous analysation was performed under the same arguments old data is removed.

        Parameters
        ----------
        parameters : List of given Arguments 

        Returns
        -------
        None.

        """
        df = pd.DataFrame( columns=['genome','width','mode','chr','biosource','tf'])
        df_length = len(df)
        df.loc[df_length] = parameters
        
        try:
            
            ori = rp().read_csv(self.path)
        
            data = ori.copy()
            
            
            data.pop('means')
            data.pop('covariances')
            data.pop('weights')
            data.pop('path')
            data.pop('time')
            data.pop('vis_filename')
            
            
            count = 0
            tochange = []
            
            print("looking for duplicates in result.csv")
            
            for row in range(len(data)):
                
                if((data.iloc[row].values == df.values).all()):
                    
                    print("replacing old item in result.csv")
                    tochange.append(count)
            
                count += 1    
            
            ori.drop(ori.index[tochange], inplace=True)
            
            os.remove(self.path)
            ori.to_csv(self.path, index=False)

        except:
            
            pass

#FOR TESTING // EXAMPLE SEE BELOW
# if __name__ == '__main__':
    
#     parameters = "hg19",	"50", "auto", "chr1", "kidney","ar"
    
#     modifyCSV(".../python-workspace/jlu-bda-2020/results/result.csv").compare(parameters)
#     #df = pd.DataFrame(data=numpy_data, columns=['genome','width','mode','chr','biosource','tf'])
