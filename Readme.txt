RLAG is the source code of the method RLAG that was proposed in the paper entitled “Identifying and ranking potential cancer drivers using representation
learning on attributed network”

Run the RLAG should under the environment of python2.7 and R >=2.7

1.invoke the main.py under the environment of python2.7 with loading the two graph of structural graph(mutated gene-gene network) and bipartite graph(mutated gene-attribute network) from the offered input data file.

2.by setting the threshold to filter out genes in each subgroup of the gene feature vectors learned in the first step, we can obtain the breast_list1 (lung_list1,  prostate_list1).Load the gene feature vectors and genes(lung,prostate,breast) in filter_subgroup.R with R.
 
3.load the list1(eg: breast_list1.Rdata),somatic mutation matrix(eg:Breast.Rdata), weight matrix of Euclidean distance(eg: breastw_dis),and genes(lung,breast,prostate) to invoke the main_function(breastw_dis)(lung,breast,prostate).
The result return a list of ranking driver genes and its corresponding score in the term of structural entropy and relative entropy.

4.use the final_scoring.R with loading the result of step3 and two data relating to genes mutations in the samples(eg:gene_sample_number and breast_sample_genelist).
The final result return a list of ranking driver genes and its corresponding score.
