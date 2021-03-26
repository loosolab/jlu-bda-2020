#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 17:40:57 2021

@author: jan
"""
from scripts.repository import Repository
from scripts.components_fit import GmFit
from scripts.visualize_data import VisualizeData as VD
from scripts.ema import EMA
from scripts.modify_csv import modifyCSV 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

class TF_analyser:
    """
    Main Class of the Analyse Part. This class is to call the whole pipe of this
    part.
    """
    
    def __init__(self, n_comps, genome, width, path, chromosome):
        """
        Initiation of Variables

        Returns
        -------
        None.

        """
        self.chr = chromosome
        self.path_results = path
            
        self.genome = genome
        self.width = width
        self.evaluate_n = True
        self.eval_size = 7
        
        if n_comps:
            
            self.evaluate_n = False
            self.n_components = n_comps
            
        path_scripts = os.path.dirname(__file__)
        path_bin = os.path.split(path_scripts)
        path_main = os.path.split(path_bin[0])
        self.path_results = os.path.join(path_main[0], 'results')
            
    # def progress(self, count, total, suffix=''):
    #     """
    #     Method to Display a progress bar
        
    #     Parameters
    #     ----------
    #     count : actual status by counted Tfs
    #     total : total Tfs
    #     suffix : TYPE, optional
    #         DESCRIPTION. The default is ''.

    #     Returns
    #     -------
    #     None.

    #     """
        
    #     bar_len = 60
    #     filled_len = int(round(bar_len * count / float(total)))
    
    #     percents = round(100.0 * count / float(total), 1)
    #     bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    #     sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    #     sys.stdout.flush()
        
    def mainloop(self, data):
        """
        Main method of the analyse part. Here the Transcription Factors 
        are analysed and a dataframe of the results is created. Also the distributions of the 
        ATAC and CHIP data is plotted and safed as png.

        Parameters
        ----------
        data: TYPE: mutiple dicts in dicts containing the actual scores in the last one
            Data to be analysed

        Returns
        -------
        resultframe: TYPE: pandas Dataframe
            result dataframe

        """
        
        
        # total = Input().number_of_chr(data)
        
        resultframe =pd.DataFrame(columns=['genome','width','mode','chr','biosource','tf','means','covariances', 'weights']) 
        
        # i = 0
        # loop all biosources
        print("analyse_main.py: unpacking")
        for biosource, b_value in data.items():
            print("unpacking: " + biosource)
            #loop all transcription factors
            for tf, tf_value in b_value.items():
                
                print('analysing: '+ tf)
                
                scoresarray = []
                #combine all scores of the chromosomes into vector-format list
                for chromosome in tf_value.values():
                
                    for array in chromosome:
    
                            scoresarray.append([array[-2], array[-1]])
            
                
                # Main().progress(i, total, '')
                
                
               # scaled_scores = TF_analyser.scale(self, scoresarray)
                distribution = np.array(scoresarray)
                
                mode = 'manual'
                
                if self.evaluate_n == True:
                    
                    mode = 'auto'
                    print("mode auto: number of components is evaluated (components_fit.py)")
                    #automated number of components evaluation  
                    all_diffs = GmFit.getDifference(self, distribution, self.eval_size)
                    self.n_components = GmFit.evaluate(self, all_diffs)
                    #plt.plot(all_diffs)
                
                single_result = EMA().emAnalyse(distribution, self.n_components)
               
                single_result.insert(0,'tf',tf)
                single_result.insert(0,'chr',self.chr)
                single_result.insert(0,'biosource',biosource)
                single_result.insert(0, 'mode', mode)
                single_result.insert(0,'width', self.width)
                single_result.insert(0,'genome', self.genome)
                
                #visualization and saving plots 
                v= VD(self.path_results, tf, self.genome, biosource, self.chr)
                path = v.displayDensityScatter(distribution, tf)
                
                v.altitudePlot(distribution, self.n_components, tf)
                z = v.contourPlot(distribution, self.n_components, tf)
            
                #Add z axis to scoresarray:
                for i in range(0,len(z)):
                    scoresarray[i].append(z[i])
                    
                #save data
                np.savetxt(path + '/' + tf + '.csv', scoresarray, delimiter=',')
                
                single_result.insert(9, 'path', path)
                
                resultframe = pd.concat([resultframe, single_result])
                print (tf + "    Done")
                # i += 1
                # Main().progress(i, total, '')
                parameters = [self.genome, self.width, mode, self.chr, biosource, tf]
                modifyCSV(self.path_results + "/result.csv").compare(parameters)
                

        #Save resultframe
        Repository().save_csv(resultframe)
            
        return resultframe
    
    def scale(self, scoresarray):
        """
        NOT USED IN THE FINAL VERSION
        =============================
        
        Method to scale the data to values from 0 to 100 

        Parameters
        ----------
        scoresarray: TYPE: list of float64 vectors
            distribution

        Returns
        -------
        scaled : TYPE: list of float64 vectors
            scaled distribution

        """
        
        count = 0

        max_x = 0
        max_y = 0
        scaled = []
        
        for i in scoresarray:
            
            if i[0] > max_x:
                max_x = i[0]
            
            if i[1] > max_y:
                max_y = i[1]
                
        for k in scoresarray:
            
            scaled_x = ((scoresarray[count][0])/max_x)*100
            scaled_y = ((scoresarray[count][1])/max_y)*100
            
            scaled.append([scaled_x, scaled_y])
            
            count  += 1
            
        return scaled
            
        
if __name__ == '__main__':

    data = Repository().inputHandler(path='/home/jan/python-workspace/angewendete_daten_analyse/testsets/calculated_data_3.pickle')
    path_scripts = os.path.dirname(__file__)
    path_bin = os.path.split(path_scripts)
    path_main = os.path.split(path_bin[0])
    path = os.path.join(path_main[0], 'results')
    
    resultframe = TF_analyser(None, "Genome", "width", path, "chr1").mainloop(data)
    # dirname = os.path.dirname(__file__)
    # resultframe.to_csv(dirname + '/result.csv', index = False, decimal=(','))
    print(resultframe)
        
                
                
            
        
