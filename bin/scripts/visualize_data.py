#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 18:15:37 2021

@author: jan
"""

import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
from sklearn.mixture import GaussianMixture
import os

class VisualizeData:
    
        
        def __init__(self,path,tf_id, genome, biosource, chromosome):
            """
            Initialize variables and set up directory if necessary

            Parameters
            ----------
            path: TYPE: str
                Path to results
            tf_id : TYPE: str
                ID of the transcription factor

            Returns
            -------
            None.

            """
            self.chromosome = chromosome
            #genome_path = (os.path.join(path, genome))
            self.path_plots = (os.path.join(path,'plots',genome ,biosource ,tf_id))
            path_scripts = os.path.dirname(__file__)
            path_bin = os.path.split(path_scripts)
            path_main = os.path.split(path_bin[0])
            self.path_visualization = os.path.join(path_main[0], "visualization","assests","img", biosource +"_"+ genome, tf_id + "_" + "".join(chromosome))
            try:
                os.makedirs(self.path_plots)
            except:
                pass
            
            try:
                os.makedirs(self.path_visualization)
            except:
                pass
            
        def makeArray(self, scores_array):
            """
            

            Parameters
            ----------
            scores_array : TYPE: list of float64
                Distribution

            Returns
            -------
            x : TYPE: nparray of float64
                Distribution
            y : TYPE: nparray of float64
                Distribution

            """
            
            x = []
            y = []
            
            for i in range(0, len(scores_array)):
                v = scores_array[i]
                xi = v[0]
                yi = v[1]
                
                x.append(xi*100)
                y.append(yi*100)
            
            np.array(x)
            np.array(y)
            
            return x,y
        
        
        #Make Density Scatter Heatmap
        def displayDensityScatter(self,scores_array, tf_id):
            """
            Method to illustrate distribution via Density Scatter(Heat-Map)

            Parameters
            ----------
            scores_array: TYPE: list of float64 vectors

            Returns
            -------
            Path: TYPE: str
                path to plots

            """
            x,y = VisualizeData.makeArray(self, scores_array)
            
            # Calculate the point density
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            
            fig, ax = plt.subplots()
            ax.scatter(x, y, c=z, s=50, edgecolors='face')
            
            ax.set(xlim=(0,100), ylim=(0,100))
            plt.xlabel("ATAC")
            plt.ylabel("Chip")
            # plt.colorbar()
            figure_path = self.path_plots + "/DensityScatter_" + tf_id + ".svg"
            plt.savefig(figure_path, format="svg")
            plt.show()
            
            return self.path_plots
        
        #Make contourPlot
        def contourPlot(self, scores_array, n_cgauss, tf_id):
            """
            Method to display distribution via contour-plot

            Parameters
            ----------
            scores_array: TYPE: list of float64 vectors
                distribution to plot
            n_cgauss: TYPE: int
                number of componets for a Gasussian Mixture Model
            Returns
            -------
            None.

            """
            
            gmm = GaussianMixture(n_components=n_cgauss)
            gmm.fit(scores_array)
            
            x,y = VisualizeData.makeArray(self, scores_array)
            # Calculate the point density
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            
            # Make the plot
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.set_ylabel('CHIP')
            ax.set_xlabel('ATAC')
            ax.plot_trisurf(x, y, z, cmap=plt.cm.coolwarm, linewidth=1, antialiased=False)
            # ax.plot_surface(x, y, z, color='b')
            figure_path = os.path.join(self.path_plots, "Contour_" + tf_id + ".svg")
            plt.savefig(figure_path, format="svg")
            vil_fig_path = os.path.join(self.path_visualization, "Contour_" + tf_id + "_" + "".join(self.chromosome) +".svg")
            plt.savefig(vil_fig_path, format="svg")
            plt.show()
            
            return z
            
        #Make altitude Plot
        def altitudePlot(self, data, n_cgauss, tf_id):
            """
            Method to display a distribution via altitude-plot. Lines for altitude 
            measurments.

            Parameters
            ----------
            data: TYPE: list of float64 vectors 
                distribution to plot
            n_cgauss: TYPE: int
                number of components 

            Returns
            -------
            None.

            """
        
            gmm = GaussianMixture(n_components=n_cgauss)
            gmm.fit(data)
            
            X, Y = np.meshgrid(np.linspace(start = -1, stop = 100, num = 100), np.linspace(start = -1, stop = 100, num = 100))
            XY = np.array([X.ravel(), Y.ravel()]).T
            Z = gmm.score_samples(XY)
            Z = Z.reshape(100,100)
    
            plt.contour(X,Y,Z)
            plt.scatter(data[:,0], data[:,1])
            
            
            figure_path = self.path_plots + "/Altitude_" + tf_id + ".svg"
            plt.savefig(figure_path, format= "svg")
            #plt.savefig(/visualization/assests/img/, kwargs)
            
            plt.show()
            
if __name__ == '__main__':
    path = "/home/python/"
    tf_id = "1234"
    genome = "genom"
    biosource = "Dito"
    
    v = VisualizeData(path, tf_id, genome, biosource, "chr1")
