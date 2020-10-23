# MultiNet
![Multinet1](Multinet1.png)

## Overview
MultiNet is a computational algorithm of topological data analysis (TDA) that manages high-dimensional datasets obtaining multiple results. It is a compound of different R functions of data analysis, having as basis the TDA's idea of representing the shape of the data cloud creating simplicial complexes from data. 

## Installation and Sample Code
Please follow the program guide in the MultiNet repository called *MultiNet_Main.R*.

## Usage 
MultiNet was designed to analyze DNA methylation arrays (450K Illumina) but it could be easily extended to other types of high-dimensional data. It provides a great quantity of results focused on decoding the correlation structure and the data patterns based on different sample groups.

It is designed to be useful in a global and local enviroment, so different analyses could be done using the addecuate MultiNet parameters, extracting multiple perspectives from the data.

## Functionalities
The *multinet* function was included in the program *MultiNet_Supplemental1.R*. 

Additionally, we have developed several auxiliary functions to program helpful graphics that allow an easy data interpretation (in this case, biological interpretation), plus output documents as spreadsheets with the results. They are collected in the program *MultiNet_Supplemental2.R* and were designed to:

1. Prepare the data, select the parameters, and generate MultiNet networks for each sample group (normally divided in overall/control/case), including colored graphs.

2. Differentiation of MultiNet networks.
    
3. Select substantial information from the networks, as DMSs, HCSs, or SMSs.
    
4. Do hierarchical cluster analysis over those DMSs and represent them in a heatmap for sample differentiation from distinct methylation patterns.
    
5. Select significantly differentiated sites over case/control groups with logistic regression. Select the most relevant CpG markers thought random forest.
    
6. Validate the significant sites found in the previous step in a different dataset (selected by the user). Select the most relevant CpG markers for sample prediction.
    
7. Study the correlation of the detected DMSs.
    
8. Analyze the regions, genes and biological pathways associated to the significant DMSs.
    
10. Apply local MultiNet over the chromosomes of interest, as those that contain a greater number of significant DMSs.
    
11. Record main results in a spreadsheet.

## Parameters
After the correct pre-processing, *multinet* requires from the following parameters to obtain the mentioned results:

1. The creation of biological results 1=Yes; 0=No

2. The creation of case: 1=Yes; 0=No

3. The creation of control: 1=Yes; 0=No

4. The initial CpG to be selected in the input array

5. The last CpG to be selected in the input array

6. The window size 

7. The window's overlap in percentage

8. The number of intervals per window

9. The interval's overlap in percentage

10. The number of cluster per interval

11. The edge-joining parameter

12. The edge-joining parameter 

13. The k-means methodology

## Other Applications
Out of the biological field, the program *MultiNet_Ibex35.R* contains the application of the *multinet* function to the Spanish stock data. It only uses a part of the algorithm, without the biological results (first parameter = 0).

## Citation
As this function and the related ones are part of my thesis work, the citation to it will be needed. Once results are published, I will include here the corresponding citation.

## Enjoy!
Select the parameteres and prepare your folder to obtain results as fast as the speed of light!
