# Searching for a new probability distribution for modeling non-scale-free heavy-tailed real-world networks

In this study, we consider large-scale network data sets from different disciplines, namely social networks, collaboration networks, web graphs, citation networks, biological networks, product co-purchasing networks, temporal networks, communication networks, ground-truth networks, and brain networks. We study several individual data sets from each discipline. These data sets are publicly available at http://snap.stanford.edu/data/index.html

Usage of the repository for the paper "Searching for a new probability distribution for modeling non-scale-free heavy-tailed real-world networks": 

1. In this repository, we present an example with ego-Twitter(In) data set (ego_twitter.csv file containing the degree fequency data). The data set can be directly imported to R statistical software for the analysis of degree distribution of networks data. 

2. The "models.R" file contains the implementation of popularly-used degree distributions, namely Lomax, power-law, power-law with cutoff, Log-normal, and Exponential distributions. Furthermore, the file "models.R" also contains the implementations of our proposed "Generalized Lomax" family of distributions, namely GLM Type-I, GLM Type-II, GLM Type-III and GLM Type-IV models. The decsirptions of all these models are provided in the manuscript titled "Searching for a new probability distribution for modeling non-scale-free heavy-tailed real-world networks". 

3. Once the implementation is done, the predicted outputs of Lomax, power-law, power-law with cutoff, Log-normal, Exponential, GLM Type-I, GLM Type-II, GLM Type-III and GLM Type-IV models are restored in a ego_twitter_output.csv file. This file is further used for the computation of different metrics for finding predictive accuracy of several models in the manuscript. 

4. Using the outputs of ego_twitter_output.csv file, we obtain the graphs (Plots of degree distributions along with different proabbility distributions) for our paper and the codes are given in figures_plots.m (MATLAB implementation file). 

5. Reults obtained in the paper for ego_twitter data can directly be computed along with the graphs and figures using the implementation files and data sets (along with outputs) given in this repository for replicability and sake of reproducibility of our paper. 
