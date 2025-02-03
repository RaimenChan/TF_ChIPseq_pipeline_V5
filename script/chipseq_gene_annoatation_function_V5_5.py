import numpy as np

def type1_downstream_geneannotation(gff, Chr, summit,distance, type3, log): 
    ###基因重叠，一个peak可能有两个type3
    if(type3!=None):
        type3 = type3.split(', ')
    else:
        type3 = []
    
    # peak 的下游没有基因， min_start ='flag'; 有基因，min_start会变成数值
    min_start = 'flag'
    
    gff_chr = gff[gff['chr']==Chr]
    gff_chr_down = gff_chr[gff_chr['start']>=summit]
    gff_chr_down = gff_chr_down[gff_chr_down['strand'] == '+']
    gff_chr_down = gff_chr_down.sort_values('start', ascending=True)
    
    ##找出不在gene body上的最近的基因
    for i in range(len(gff_chr_down)):
        tem_gene = gff_chr_down.iloc[i]['ID']
        tem_gene = tem_gene.split(';')[0]
        tem_gene = tem_gene.split('-')[0]
        tem_gene = tem_gene.split('=')[1]
        if(tem_gene in type3):
            continue
        else:
            min_start = gff_chr_down.iloc[i]['start']
            break
    
    #peak的下游没有基因
    if(min_start=='flag'):
        return None, None
    
    d = min_start - summit
    if(d<=distance):
        gff_chr_down_target = gff_chr_down[gff_chr_down['start']==min_start]
        if(len(gff_chr_down_target)==1):
            target_ID_value = gff_chr_down_target.iloc[0]['ID']
            target_ID = target_ID_value.split(';')[0]
            target_ID = target_ID.split('-')[0]
            target_ID = target_ID.split('=')[1]
            return target_ID, d
        else:
            print(Chr,'\t', summit, 'down stream abnormal:', "type1 down more than 1 genes")
            log.write(Chr + " " + str(summit) + " " + 'down stream abnormal: ' + "type1 down more than 1 genes\n")
            for i in range(len(gff_chr_down_target)):
                print(gff_chr_down_target.iloc[i]['ID'])
                log.write(gff_chr_down_target.iloc[i]['ID'])
                log.write("\n")
            print()
            log.write("\n")
            



            target_IDs = []
            for i in range(len(gff_chr_down_target)):
                target_ID_value = gff_chr_down_target.iloc[i]['ID']
                target_ID = target_ID_value.split(';')[0]
                target_ID = target_ID.split('-')[0]
                target_ID = target_ID.split('=')[1]
                target_IDs.append(target_ID)
            
            gene_IDs = ''
            for i,ID in enumerate(target_IDs):
                gene_IDs = gene_IDs + ID
                if(i!=len(target_IDs)-1):
                    gene_IDs = gene_IDs + ', '
            return gene_IDs, d
    else:
        return None, None

    
def type1_upstream_geneannotation(gff, Chr, summit,distance, type3, log):
    if(type3 != None):
        type3 = type3.split(', ')
    else:
        type3 = []
    
    # peak 的上游没有基因, max_end = 'flag'; 有基因，max_end 会变成数值
    max_end = 'flag'
    
    gff_chr = gff[gff['chr']==Chr]
    gff_chr_up = gff_chr[gff_chr['end']<=summit]
    gff_chr_up = gff_chr_up[gff_chr_up['strand'] == '-']
    gff_chr_up = gff_chr_up.sort_values('end', ascending=False)
    
    ##找出不在gene body上的最近的基因
    for i in range(len(gff_chr_up)):
        tem_gene = gff_chr_up.iloc[i]['ID']
        tem_gene = tem_gene.split(';')[0]
        tem_gene = tem_gene.split('-')[0]
        tem_gene = tem_gene.split('=')[1]
        if(tem_gene in type3):
            continue
        else:
            max_end = gff_chr_up.iloc[i]['end']
            break
    
    # peak 的上游没有基因
    if(max_end == 'flag'):
        return None, None
    
    d = summit - max_end
    if(d<=distance):
        gff_chr_up_target = gff_chr_up[gff_chr_up['end']==max_end]
        if(len(gff_chr_up_target)==1):
            target_ID_value = gff_chr_up_target.iloc[0]['ID']
            target_ID = target_ID_value.split(';')[0]
            target_ID = target_ID.split('-')[0]
            target_ID = target_ID.split('=')[1]
            return target_ID, d
        else:
            print(Chr,'\t', summit, 'up stream abnormal:', "type1 up more than 1 genes")
            log.write(Chr + " " + str(summit) + " " +  'up stream abnormal: ' + "type1 up more than 1 genes\n")
            for i in range(len(gff_chr_up_target)):
                print(gff_chr_up_target.iloc[i]['ID'])
                log.write(gff_chr_up_target.iloc[i]['ID'])
                log.write("\n")
            print()
            log.write("\n")



            target_IDs = []
            for i in range(len(gff_chr_up_target)):
                target_ID_value = gff_chr_up_target.iloc[i]['ID']
                target_ID = target_ID_value.split(';')[0]
                target_ID = target_ID.split('-')[0]
                target_ID = target_ID.split('=')[1]
                target_IDs.append(target_ID)
            
            gene_IDs = ''
            for i,ID in enumerate(target_IDs):
                gene_IDs = gene_IDs + ID
                if(i!=len(target_IDs)-1):
                    gene_IDs = gene_IDs + ', '
            return gene_IDs, d
    else:
        return None, None


# binding in 3'-UTR or 5'-UTR
def type2_geneannotation(gff, Chr, summit):
    gff_chr = gff[gff['chr']==Chr]
    gff_chr = gff_chr[gff_chr['start']<=summit]
    gff_chr = gff_chr[gff_chr['end']>=summit]
    if(len(gff_chr)==1):
        target_ID_value = gff_chr.iloc[0]['ID']
        target_ID = target_ID_value.split('-')[0]
        target_ID = target_ID.split('=')[1]
        target_ID  = target_ID .split('_')[1]
        return target_ID
    
    elif(len(gff_chr)==0):
        return None
    
    #
    else:
        print(Chr,' ', summit, ' ', "type2 more than 1 gene(five_prime_UTR)")
        for i in range(len(gff_chr)):
            print(gff_chr.iloc[i]['ID'])
        print()
        
        target_IDs = []
        for i in range(len(gff_chr)):
            target_ID_value = gff_chr.iloc[i]['ID']
            target_ID = target_ID_value.split('-')[0]
            target_ID = target_ID.split('=')[1]
            target_ID  = target_ID .split('_')[1]
            target_IDs.append(target_ID)
        
        gene_IDs = ''
        for i,ID in enumerate(target_IDs):
            gene_IDs = gene_IDs + ID
            if(i!=len(target_IDs)-1):
                gene_IDs = gene_IDs + ', '
        return gene_IDs



# this function will return genes also in type2, need to remove the type2  gene.
def type3_geneannotation(gff, Chr, summit, log):
    gff_chr = gff[gff['chr']==Chr]
    gff_chr = gff_chr[gff_chr['start']<=summit]
    gff_chr = gff_chr[gff_chr['end']>=summit]
    
    # peak 不绑定在gene body上
    if(len(gff_chr)==0):
        return None
    
    # peak 绑定在gene body上，正常情况
    elif(len(gff_chr)==1):
        target_ID_value = gff_chr.iloc[0]['ID']
        target_ID = target_ID_value.split(';')[0]
        target_ID = target_ID.split('=')[1]
        return target_ID
    
    # peak绑定在gene body上，但大于一个基因，异常情况
    else:
        print(Chr,' ', summit, ' type3 more than 1 gene')
        log.write(Chr + ' ' + str(summit) + " " + 'type3 more than 1 gene\n')

        for i in range(len(gff_chr)):
            print(gff_chr.iloc[i]['ID'])
            log.write(gff_chr.iloc[i]['ID'])
            log.write("\n")
        print()
        log.write("\n")


        target_IDs = []
        for i in range(len(gff_chr)):
            target_ID_value = gff_chr.iloc[i]['ID']
            target_ID = target_ID_value.split(';')[0]
            target_ID = target_ID.split('=')[1]
            target_IDs.append(target_ID)
        
        gene_IDs = ''
        for i,ID in enumerate(target_IDs):
            gene_IDs = gene_IDs + ID
            if(i!=len(target_IDs)-1):
                gene_IDs = gene_IDs + ', '
        return gene_IDs


def remove_type3_in_type2(three, five, gene_body):
    try:
        gene_body = gene_body.split(", ")
    except:
        gene_body = []
    try:
        three = three.split(", ")
    except:
        three = []
    try:
        five = five.split(", ")
    except:
        five = []
    out = ""
    for gene in gene_body:
        if(gene in three):
            continue
        elif(gene in five):
            continue
        else:
            if(out == ""):
                out = gene
            else:
                out = out + ", " + gene

    if(out == ""):
        return None
    else:
        return out
    
    

def target_gene_list(type3, type1_down, type1_up):
    target_genes =[]
    for genes in type3:
        target_genes = target_genes + genes.split(',')
    for genes in type1_down:
        target_genes = target_genes + genes.split(',')
    for genes in type1_up:
        target_genes = target_genes + genes.split(',')
    return target_genes


def two_type1_selection(d_ratio, downg, down_d, upg, up_d, d2):
    ##downg:down stream type1 gene
    ##upg: upstream type1 gene
    ##down_d, up_d: up/down distance to peak summit

    #只有上下游都有基因的时候再根据d2筛选,不然输入什么，输出什么
    if(downg != None and upg != None):
        
        ## 上下游都有基因，而且两个基因TSS/ATG之间的距离小于d2
        if(down_d + up_d <=d2):
            return downg, down_d, upg, up_d

        #peak 靠近上游的基因
        if(down_d >= up_d):
            ratio = down_d / up_d
            if(ratio < d_ratio):
                return downg, down_d, upg, up_d
            else:
                return None, float('nan'), upg, up_d
        
        #peak 靠近下游的基因  
        if(down_d < up_d):
            ratio = up_d / down_d
            if(ratio < d_ratio):
                return downg, down_d, upg, up_d
            else:
                return downg, down_d, None, float('nan')
            
    else:
        return downg, down_d, upg, up_d
    

def summit_signal(bdg_file, Chr, summit):
    
    tem_bdg = bdg_file[bdg_file['chr']==Chr]
    
    tem_bdg2 = tem_bdg[ (tem_bdg['start'] <= summit) & (tem_bdg['end'] >= summit)]
    
    signal = tem_bdg2.iloc[0,3]

    return signal


def add_description(ID, ID_description):
    if(ID==None):
        return '#N/A'
    else:
        IDs = ID.split(', ')
        descriptions =''
        for id in IDs:
            description = ID_description[ID_description['ID']==id]['description'].values[0]
            descriptions = descriptions + description +", "
        return descriptions
    



# heatmap function get the binding signal
def get_bdg_200bp_per_bins_promoter(bdg, Chr, start, end):
    tem_bdg = bdg[bdg['chr']==Chr].copy()
    
    tem_bdg_start = tem_bdg[tem_bdg['start']<start]
    start_threshold = max(tem_bdg_start['start'])
    
    tem_bdg_end = tem_bdg[tem_bdg['end']>=end]
    end_threshold = min(tem_bdg_end['end'])
    
    tem_bdg = tem_bdg[tem_bdg['start']>=start_threshold]
    tem_bdg = tem_bdg[tem_bdg['end']<=end_threshold]
    
    tem_bdg.iloc[0,1]=start
    tem_bdg.iloc[-1,2]=end
    #下一行的start和前一行的end为同一个数，取前一行的end
    scores = []
    for i in range(len(tem_bdg)):
        length = tem_bdg.iloc[i,2] - tem_bdg.iloc[i,1]
        score = tem_bdg.iloc[i,3]
        List = [score] * length

        scores = scores + List
    merge_scores=[]
    for i in range(0,2000,100):
        merge_scores.append(np.mean(scores[i:i+100]))
    
    return merge_scores

def get_bdg_200bp_per_bins_tail(bdg, Chr, start, end):
    tem_bdg = bdg[bdg['chr']==Chr].copy()
    
    tem_bdg_start = tem_bdg[tem_bdg['start']<start]
    start_threshold = max(tem_bdg_start['start'])
    
    tem_bdg_end = tem_bdg[tem_bdg['end']>=end]
    end_threshold = min(tem_bdg_end['end'])
    
    tem_bdg = tem_bdg[tem_bdg['start']>=start_threshold]
    tem_bdg = tem_bdg[tem_bdg['end']<=end_threshold]
    
    tem_bdg.iloc[0,1]=start
    tem_bdg.iloc[-1,2]=end
    #下一行的start和前一行的end为同一个数，取前一行的end
    scores = []
    for i in range(len(tem_bdg)):
        length = tem_bdg.iloc[i,2] - tem_bdg.iloc[i,1]
        score = tem_bdg.iloc[i,3]
        List = [score] * length

        scores = scores + List
    merge_scores=[]
    for i in range(0,2000,100):
        merge_scores.append(np.mean(scores[i:i+100]))
    
    return merge_scores

def get_bdg_gene_body(bdg, Chr, start, end):
    tem_bdg = bdg[bdg['chr']==Chr].copy()
    
    tem_bdg_start = tem_bdg[tem_bdg['start']<start]
    start_threshold = max(tem_bdg_start['start'])
    
    tem_bdg_end = tem_bdg[tem_bdg['end']>=end]
    end_threshold = min(tem_bdg_end['end'])
    
    tem_bdg = tem_bdg[tem_bdg['start']>=start_threshold]
    tem_bdg = tem_bdg[tem_bdg['end']<=end_threshold]
    
    tem_bdg.iloc[0,1]=start
    tem_bdg.iloc[-1,2]=end
    #下一行的start和前一行的end为同一个数，取前一行的end
    scores = []
    for i in range(len(tem_bdg)):
        length = tem_bdg.iloc[i,2] - tem_bdg.iloc[i,1]
        score = tem_bdg.iloc[i,3]
        List = [score] * length

        scores = scores + List
    
    merge_scores=[]
    length = end - start
    step = int(length/20)
    merge_scores=[]
    k=0
    for i in range(0,length,step):
        k=k+1
        Mean = np.mean(scores[i:i+step])
        if(k==20):
            Mean = np.mean(scores[i:])
            merge_scores.append(Mean)
            break
        merge_scores.append(Mean)
        

    return merge_scores

def extract_common_name(description):
    if("Name=" in description):
        return description.split("Name=")[1].split(";")[0]
    else:
        return "-"