# MultiNet

MultiNet is a computational algorithm of topological data analysis (TDA) that manages high-dimensional datasets obtaining multiple results. It is a compound of different R functions of data analysis, having as basis the TDA's idea of representing the shape of the data cloud creating simplicial complexes from data. 

MultiNet was designed to analyze DNA methylation arrays (450K Illumina) but it could be easily extended to other types of high-dimensional data. It provides a great quantity of results focused on decoding the correlation structure and the data patterns based on different sample groups.

It is designed to be useful in a global and local enviroment, so different analyses could be done using the addecuate MultiNet parameters, extracting multiple perspectives from the data.

MultiNet needs the input of the main dataset plus the filter functions, and it only asks for a group of parameters before creating the networks. In addition, we have developed several functions to program helpful graphics that allow an easy data interpretation (in this case, biological interpretation), plus output documents as spreadsheets with the results. They were designed to:

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

After the correct pre-processing, MultiNet requires from the following parameters to obtain the mentioned results:

1. The creation of case: 1=Yes; 0=No

2. The creation of control: 1=Yes; 0=No

3. The initial CpG to be selected in the input array

4. The last CpG to be selected in the input array

5. The window size 

6. The window's overlap in percentage

7. The number of intervals per window

8. The interval's overlap in percentage

9. The number of cluster per interval

10. The edge-joining parameter

11. The edge-joining parameter 

12. The k-means methodology


Select the parameteres and prepare your folder to obtain results as fast as the speed of light!
