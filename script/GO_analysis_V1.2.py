import pandas as pd
from scipy.stats import hypergeom
from goatools import obo_parser
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import seaborn as sns
import os
import argparse

##################################### arguments #####################################################

parser = argparse.ArgumentParser(description="According to the hypergeometric distribution, calculate the enrichment probability of the target genes in various GO terms,\n\
 and select the GO terms with a p-value less than 0.05 for visualization.\n",
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-gaf", help="gaf file or file contain the Gene2GO relation. If the file not in .gaf format, please use -gene_index\
-go_index, -aspect_index to assign the column indices of gene, go and aspect"
                                        )
parser.add_argument("-gene_index", help="gene column index in gaf file, default=1", default=1)
parser.add_argument("-GO_index", help="GO column index in gaf file, default=4", default=4)
parser.add_argument("-aspect_index", help="gene column index in gaf file, default=8", default=8)

parser.add_argument("-target_gene_file", help="target gene should save in a txt file format like \ngene1\ngene2\ngene3\n...")
parser.add_argument("-obo", help=".obo file")
parser.add_argument("-o", help="output", default='./output')

args = parser.parse_args()
gene_index = int(args.gene_index)
GO_index = int(args.GO_index)
aspect_index = int(args.aspect_index)



#################################### read input files ###############################################
gaf =pd.read_csv(args.gaf, sep='\t', comment='!',header=None)
gaf = gaf[[gene_index, GO_index,aspect_index]]
Aspects = list(set(gaf[aspect_index]))

go = obo_parser.GODag(args.obo)

with open(args.target_gene_file, 'r') as file:
    file_content = file.read()
Target_gene_list = file_content.split('\n')


#################################### enrichment analysis ##############################################
results = []
# calculate BP, CC, MF seperately
for aspect in Aspects:
    
    #GO2gene association dict
    gaf_aspect = gaf[gaf[aspect_index]==aspect].copy()
    assoc={}
    group = gaf_aspect.groupby([GO_index])
    for key,tem_df in group:
        GO_ID = key[0]
        Gene_IDs = tem_df[gene_index].to_numpy().tolist()
        Gene_IDs = list(set(Gene_IDs))
        assoc[GO_ID] = Gene_IDs
        
        
    background_gene_ID =[]
    for GO in assoc.keys():
        background_gene_ID = background_gene_ID + assoc[GO]
    
    background_gene_ID = list(set(background_gene_ID))
    
    #Target gene list 和 background 的交集
    Study = list( 
        set(Target_gene_list).intersection( set(background_gene_ID) )
    )
    
    
    
    for GO in assoc.keys():
        GO_result ={}
        Genes_in_GO = assoc[GO]
        Population = len(background_gene_ID)
        Term_size = len(Genes_in_GO)
        Intersection = list(
        set(Genes_in_GO).intersection( set(Study))
        )
        Study_size = len(Study)
        Intersection_size = len(Intersection)
        
        #Intersection_size - 1才是我们想要的pvalue
        p_value = hypergeom.sf(Intersection_size - 1, Population, Term_size, Study_size)
        
        GO_result['GO Term'] = GO
        GO_result['Category'] = aspect
        GO_result['Population'] = Population
        GO_result['Term_size'] = Term_size
        GO_result['Study_size'] = Study_size
        GO_result['Intersection_size'] = Intersection_size
        GO_result['p_value'] = p_value
        GO_result['GO_item'] = Genes_in_GO
        GO_result['Study_item'] = Study
        GO_result['Intersection'] = Intersection
        
        results.append(GO_result)
    
results_df = pd.DataFrame(results)
results_df = results_df.sort_values(['Category', 'p_value'])
## add GO info (Term name and depth to result)
Term_names = []
Term_depths = []
for i in range(len(results_df)):
    GO_term = results_df.iloc[i]['GO Term']
    if GO_term in go:
        GO_term_info = go[GO_term]
        Term_name = GO_term_info.name
        Term_depth = GO_term_info.depth
        Term_names.append(Term_name)
        Term_depths.append(Term_depth)
    else:
        Term_names.append('GO term not in obo file')
        Term_depths.append('GO term not in obo file')
results_df.insert(1,'Term_name', Term_names)
results_df.insert(2,'Term_depth', Term_depths)

results_df.to_excel('%s.xlsx'%args.o, index=False)



################################ visualization #################################
results_df['Gene%']=results_df['Intersection_size']/results_df['Study_size']*100
def Replace_No_Term_Name(Term_name, GO_term):
    if(Term_name == 'GO term not in obo file'):
        return GO_term
    else:
        return Term_name

def rename_duplicate(GO_term, Term_name, duplicate):
    if(duplicate==False):
        return Term_name
    else:
        return GO_term + '_' + Term_name
    
    
results_df['Term_name'] = results_df.apply(lambda X:Replace_No_Term_Name(X['Term_name'], X['GO Term']),axis=1)

duplicates = results_df['Term_name'].duplicated(keep = False)
results_df['duplicate'] = duplicates
results_df['Term_name'] = results_df.apply(lambda X:rename_duplicate(X['GO Term'], X['Term_name'], X['duplicate']), axis=1)
for category in Aspects:
# for category in [ 'C']:
    df_cate = results_df[results_df['Category']==category].copy()
    df_cate = df_cate.sort_values(['p_value','Term_depth','GO Term'],ascending=[True,True,True])
    df_cate_005 = df_cate[df_cate['p_value']<0.05].copy()
    
    
    ##### if too many items in each category whose p < 0.05, show them in different picture, each picture show 20 items
    for i in range(len(df_cate_005)//20+1):
        if(20*i+20<=len(df_cate_005)):
            df_cate_005_i = df_cate_005.iloc[20*i:20*i+20,:]
        else:
            df_cate_005_i = df_cate_005.iloc[20*i:len(df_cate_005),:]
        
        if(len(df_cate_005_i)==0):
            continue
        
        
        #### parameter used to draw the picture
        X = df_cate_005_i['Gene%'].to_numpy().tolist()
        Y = df_cate_005_i['Term_name'].to_numpy().tolist()
        P = df_cate_005_i['p_value'].to_numpy().tolist()
        N = df_cate_005_i['Intersection_size'].to_numpy().tolist()
        M = df_cate_005_i['Term_size'].to_numpy().tolist()
        
#         X = X[::-1]
#         Y = Y[::-1]
#         P = P[::-1]
#         N = N[::-1]
        
        
        
        
        #### draw the picture
        plt.figure(figsize=(10,17))
        
        norm = plt.Normalize(0, 0.1)
        cmap = cm.get_cmap('Blues_r')  # Replace 'cool' with the desired colormap
        colors = cmap(norm(P)).tolist()
        
        sns.barplot(x=X, y=Y, palette=colors, dodge=False)
        #plt.barh( Y, X, color=colors, height =0.8, alpha =1 )
        cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', ax=plt.gca())
        cbar.ax.tick_params(labelsize=20)  # Increase the font size of the colorbar tick labels
        cbar.set_label('p value', fontsize=20, fontweight='bold')
        
        plt.xlabel('Gene[%]', fontsize=25)
        plt.xlim(0,max(X)+10)
        plt.yticks(fontsize=20, fontname='Times New Roman', fontweight='bold')
        plt.xticks(fontsize=20, fontname='Times New Roman')
        
        for i2, (x, y, p, n, m) in enumerate(zip(X, Y, P, N, M)):
            plt.text(x+0.3, i2, 'n=%d, m=%d \np=%.2e'%(n, m, p), ha='left', va='center', font='Times New Roman', fontsize = 15)
            
        plt.savefig('%s'%args.o + '_%s'%category + '_%d.jpeg'%i, dpi=700, bbox_inches='tight')
        plt.savefig('%s'%args.o + '_%s'%category + '_%d.pdf'%i, dpi=700, bbox_inches='tight')
        #plt.show()
    

print("end")