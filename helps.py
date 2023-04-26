import numpy as np
import pandas as pd
import scipy.stats as stats
import json
import scipy
import nmi
import collections
from numpy import nan
from statsmodels.stats.multitest import fdrcorrection
from time import localtime, strftime



def save_levels_nSBM(model, dataset):
    
    	model.save_graph(filename=f"{dataset}-graph.xml.gz")  
    	model.dump_model(filename=f"{dataset}-model.pkl")
    	with open(f'{dataset}-entropy.txt', 'w') as f:
        	f.write(str(model.state.entropy()))
        
    	for l in range(len(model.state.get_levels()) - 2)[::-1]:
            print(f"------Begin level {l}------",strftime("%Y-%m-%d %H:%M:%S", localtime()))
            data=model.get_groups(l)

            #P(main_feature|main_topic)
            p_w_tw=pd.DataFrame(data=data["p_w_tw"], index=model.words, 
        	     columns=[f"main_feature_topic_{i}" for i in range(data["p_w_tw"].shape[1])])
            p_w_tw.to_csv(f"{dataset}-level-{l}-main_feature-topics.csv")
        
            #P(meta_feature|meta_topic_feature) for each type of meta_feature
            for k in range(model.nbranches):
                feat_topic=pd.DataFrame(data=data["p_w_key_tk"][k], index=model.keywords[k],
        	                   columns=[f"meta_features_{k}_topic_{i}" for i in range(data["p_w_key_tk"][k].shape[1])])
                feat_topic.to_csv(f"{dataset}-level-{l}-meta_features_{k}-topics.csv")
        
        	#P(documet|cluster)
            pd.DataFrame(data=data["p_td_d"], columns=model.documents).to_csv(f"{dataset}-level-{l}-clusters.csv")
        
            #P(main_topic|documets)
            p_tw_d=pd.DataFrame(data=data["p_tw_d"].T, index=model.documents, 
        	     columns=[f"main_feature_topic_{i}" for i in range(data["p_w_tw"].shape[1])])
            p_tw_d.to_csv(f"{dataset}-level-{l}-main_topics-documents.csv")
        	
            #P(meta_topic|document) for each kind of meta_feature
            for k in range(model.nbranches):
                p_tk_d=pd.DataFrame(data=data["p_tk_d"][k].T, index=model.documents, 
                	columns=[f"meta_features_{k}_topic_{i}" for i in range(data["p_w_key_tk"][k].shape[1])])
                p_tk_d.to_csv(f"{dataset}-level-{l}-meta_topics_{k}-topics-documets.csv")      
        	


def hSBM_data(dataset, labels, info, lev):
    '''
    It takes the files from hSBM and create more readable files
    
    dataset (string) = experiment name
    labels (dataframe) = indexed by the samples' names
    info (dataframe) = indexed by the ensamble code of each gene and with one columns called "Gene name". It's used to translate the ENSG code to the gene name.
    level (int) = level of hSBM you want to analyse
    '''
    
    outcome={}
    path=f"Results/{dataset}/hSBM"
    
    #P(sample|cluster)
    with open(f"{path}/{dataset}-cluster-level-{lev}.txt") as f:
        clusters=json.load(f)
    cluster_df=pd.DataFrame.from_dict(clusters,orient="index")
    labels["hSBM"]="--"
    for i in range(len(clusters)):
        labels["hSBM"].loc[np.asarray(np.asarray(clusters[str(i)])[:,0])]=i
    labels.hSBM=labels.hSBM.astype(str)
 
    #P(topic|sample)
    p_topic_sample=pd.read_csv(f"{path}/{dataset}-topsbm_level_{lev}_topic_dist.csv",index_col=0)
    p_topic_sample.index=p_topic_sample.index.astype(str)
    p_topic_sample.columns=[(f"Topic {i+1}") for i in range(len(p_topic_sample.columns))]
    p_topic_sample.index.rename("Sample",inplace=True)
    
    #Pc(topic|sample)
    p_topic_sample_cent=p_topic_sample-p_topic_sample.mean(axis=0)

    #P(gene|topic)
    with open(f"{path}/{dataset}-topics-level-{lev}.txt") as f:
        topics=json.load(f)
    topics=pd.DataFrame.from_dict(topics,orient="index")
    topics_name=[(f"Topic {i+1}") for i in range(len(topics))]

    temp=np.zeros((len(topics),len(topics.columns)))
    temp1=[]
    for i in range(len(topics)):
        temp1.append([])
        for j in range(len(topics.iloc[i].dropna())):
            temp[i][j]=topics.iloc[i][j][1]
            temp1[i].append(topics.iloc[i][j][0])

    topic_gene_prob=pd.DataFrame(temp.T, columns=topics_name)
    topic_gene=pd.DataFrame(temp1, index=topics_name).T
    topic_gene_prob.index.rename("P(gene|topic)",inplace=True)
    topic_gene_prob.replace(to_replace=0,value=nan,inplace=True)
    
    topic_gene_genename=pd.DataFrame(index=topic_gene.index, columns=topic_gene.columns)
    for topic in topic_gene_genename.columns:
        topic_gene_genename[topic]=ens_to_name_l(topic_gene[topic], info)
    	
    topic_gene_raw=pd.DataFrame(index=flat_list([topic_gene_genename[top].dropna() for top in topic_gene_genename.columns]), columns=topic_gene.columns)
    for top in topic_gene_raw.columns:
        topic_gene_raw[top].loc[list(topic_gene_genename[top].dropna())]=list(topic_gene_prob[top].dropna())
        
    
    #Return 7 dataframes:
    outcome["sample_cluster"]=labels
    outcome["topic_gene_prob"]=topic_gene_prob
    outcome["topic_gene"]=topic_gene
    outcome["topic_gene_raw"]=topic_gene_raw
    outcome["topic_gene_genename"]=topic_gene_genename
    outcome["p_topic_sample_cent"]=p_topic_sample_cent
    outcome["p_topic_sample"]=p_topic_sample
    return outcome
                    
                  
def topic_gene_nSBM(file, top_name, info, level, dataset):
    '''
    It transforms the raw file with topic-gene probability distribution into four useful files
    
    file (string) = path of the file containing the topic_gene association
    top_name (string) = name of the topic family (e.g. "mRNA-Topic", "lncRNA-Topic", etc.)
    info (dataframe) = indexed by the ensamble code of each gene and with one columns called "Gene name". It's used to translate the ENSG code to the gene name.
    level (int) = level of nSBM you want to analyse
    dataset (string) = experiment name
    '''

    topic_gene=pd.read_csv(file, index_col=0)
    if topic_gene.index[0][0]=="#":
        topic_gene.index=[ind[1:] for ind in topic_gene.index]
        
    topic_gene.replace(to_replace=0, value=np.nan, inplace=True)
    topic_gene.columns=[f"{top_name}-Topic {i}" for i in range(len(topic_gene.columns))]

    topic_gene_raw=pd.DataFrame(index=topic_gene.index, columns=topic_gene.columns, data=topic_gene.values)
    topic_gene_raw.index=ens_to_name_l(topic_gene_raw.index, info)
    
    tgp={}
    tg={}
    for feat in topic_gene.columns:
        tgp[feat]=sorted(np.array(topic_gene[feat].dropna()), reverse=True)
        tg[feat]=topic_gene[feat].dropna().sort_values(ascending=False).index
    topic_gene=pd.DataFrame.from_dict(tg,orient="index").T
    topic_gene_prob=pd.DataFrame.from_dict(tgp,orient="index").T
    topic_gene.columns=[f"{top_name}-Topic {i}" for i in range(len(topic_gene.columns))]
    topic_gene.to_csv(f"Results/{dataset}/nSBM/Data/{dataset}-level-{level}-{top_name}-topic-gene.csv")

    topic_gene_prob.columns=topic_gene.columns
    topic_gene_prob.to_csv(f"Results/{dataset}/nSBM/Data/{dataset}-level-{level}-{top_name}-topic-gene-prob.csv")
    
    topic_gene_genename=pd.DataFrame(index=topic_gene.index, columns=topic_gene.columns, data=topic_gene.values)
    for col in topic_gene_genename.columns:
    	topic_gene_genename[col]=ens_to_name_l(topic_gene_genename[col], info)
    topic_gene_genename.to_csv(f"Results/{dataset}/nSBM/Data/{dataset}-level-{level}-{top_name}-topic-gene-genename.csv")
    
    return topic_gene, topic_gene_prob, topic_gene_genename, topic_gene_raw


def p_c_t_c_nSBM(file, top_name, labels, level, dataset):
    '''
    It takes the p_topic_cell distribution and builds the p_centered_topic_cell and the p_centered_topic_group for an arbitrary partition of cells
    
    file (string) = path of the file containing the topic_cell association
    top_name (string) = name of the topic family (e.g. "mRNA-Topic", "lncRNA-Topic", etc.)
    params (list) = the first element contains the list with the names of the partition (e.g. ["cluster 1", "cluster 2", ...]) and the second element the "name" of the group that has
                    to correspond to one column in the dataframe labels
    labels )dataframe) = indexed by the samples' names. It must contain a column whose name correspon to params[1]
    dataset (string) = experiment name
    level (int) = level of nSBM you want to analyse
    '''
    clusters=sorted(list(set(labels.nSBM)))
    
    p_topic_cell=pd.read_csv(file, index_col=0)
    p_topic_cell.columns=[f"{top_name}-Topic {i}" for i in range(len(p_topic_cell.columns))]

    p_topic_cell.to_csv(f"Results/{dataset}/nSBM/Data/{dataset}-level-{level}-{top_name}-p_topic_sample.csv")
    p_c_topic_cell=p_topic_cell-p_topic_cell.mean()

    p_c_topic_class=pd.DataFrame(index=clusters, columns=p_c_topic_cell.columns)
    for cla in clusters:
        insample=labels[labels["nSBM"]==cla].index
        p_c_topic_class.loc[cla]=p_c_topic_cell.loc[insample].mean()
    p_c_topic_class
    return p_c_topic_class, p_c_topic_cell, p_topic_cell
    
    
def topic_up(group, p_t_c_c, n_sigma):
    '''
    Function to select the enriched topic for each group of cells.
    
    group (list) = contains the names of the possibile group a cell may belongs (e.g. ["cluster 1", "cluster 2", ...] 
    p_t_c_c (dataframe) = dataframe indexed by the group and the columns to correspond to the the topics' names. Each entry is a float and shows the P_cenetred(topic|group)
    n_sigma (float) = threshold used to exclude topics
    '''
    
    if group != p_t_c_c.index.to_list():
        print("Class != index")
        return None
    else:
        topics={}
        for cla in group:
            m=np.mean(p_t_c_c.loc[cla].mean())
            s=np.mean(p_t_c_c.loc[cla].std())
            tops=[]
            for topic in p_t_c_c.columns:
                if p_t_c_c.loc[cla][topic] >= m + n_sigma * s:
                    tops.append(topic)
            topics[cla]=tops
    return topics


def topic_dn(group, p_t_c_c, n_sigma):
    '''
    Function to select the repressed topic for each group of cells.
    
    group (list) = contains the names of the possibile group a cell may belongs (e.g. ["cluster 1", "cluster 2", ...] 
    p_t_c_c (dataframe) = dataframe indexed by the group and the columns to correspond to the the topics' names. Each entry is a float and shows the P_cenetred(topic|group)
    n_sigma (float) = threshold used to exclude topics
    '''
    
    if group != p_t_c_c.index.to_list():
        print("Class != index")
        return None
    else:
        topics={}
        for cla in group:
            m=p_t_c_c.loc[cla].mean()
            s=p_t_c_c.loc[cla].std()
            tops=[]
            for topic in p_t_c_c.columns:
                if p_t_c_c.loc[cla][topic] <= m - n_sigma * s:
                    tops.append(topic)
            topics[cla]=tops
    return topics


def check_diz(diz):
    '''
    It checks if a dictionary has a empty key of or two values belonging to two different keys are equal
    '''
    for key in diz.keys():
        for k in diz.keys():
            if key==k:
                continue
            else:
                if diz[key]==diz[k] or len(diz[key])==0:
                    return True
                else: 
                    continue
    return False


def loop_topics(group, p_t_c_c, direction="up"):
    '''
    Loop over all the possibile thresholds to select enriched/repressed topics checking each step if the stop condition is fullfilled. 
    
    group (list) = contains the names of the possibile group a cell may belongs (e.g. ["cluster 1", "cluster 2", ...])
    p_t_c_c (dataframe) = dataframe indexed by the group and the columns to correspond to the the topics' names. Each entry is a float and shows the P_centered(topic|group)
    direction (string) = define if you demand for enriched (up) or repressed (down) topics
        
    '''

    if direction=="up":
        for i in np.linspace(3,0,61):
            diz = topic_up(group, p_t_c_c, n_sigma=i)
            if check_diz(diz)==False:
                print("n_sigma=",i)
                return diz, i
            else: 
                continue
    
    if direction=="down":

        for i in np.linspace(3,0,61):
            diz = topic_dn(group, p_t_c_c, n_sigma=i)
            if check_diz(diz)==False:
                print("n_sigma=",i)
                return diz, i
            else: 
                continue    
 
    
def enrichment_test(database, topics, p_g_t, info, top_name, path_to_save, dataset, level):
    '''
    It performs the enrichement test (through an hypergeometric test) for a list of topics given a specific database.
    
    database (dictionary) = each key correspond to the name of the gene set and each value is a list of genes
    p_g_t (dataframe) = P(gene|topic)
    info (dataframe) = indexed by the ensamble code of each gene and with one columns called "Gene name". It's used to translate the ENSG code to the gene name.
    top_name (string) = name of the topic family (e.g. "mRNA-Topic", "lncRNA-Topic", etc.)
    path_to_save (string) = path where you desire to save the results
    dataset (string) = experiment name
    level (int) = level of hSBM you want to analyse  
    
    '''

    M=len(flat_list([database[k] for k in database.keys()]))
    print(len(database.keys()), M)
    hyperg_val={}
    hyperg_name={}
    inter_genes={}
    len_inter_genes={}
    p_gene_topic={}
    for topic in topics.columns:
        print(topic)
        hyperg_val[topic]=[]
        hyperg_name[topic]=[]
        inter_genes[topic]=[]
        len_inter_genes[topic]=[]
        p_gene_topic[topic]=[]
        my_list=ens_to_name_l(topics[topic].dropna(), info)
        for l in database.keys():
            hyperg_val[topic].append(hypergeom(my_list,database[l], M=M))
            hyperg_name[topic].append(my_list)
            inter_genes[topic].append(intersection([my_list,database[l]]))
            len_inter_genes[topic].append(len(intersection([my_list,database[l]])))
            p_gene_topic[topic].append(list(p_g_t.loc[intersection([my_list,database[l]])][topic]))
        result=pd.DataFrame(index=database.keys(),data=hyperg_val[topic],columns=["Pvalue"])
        result["fdr"]=-np.log10(fdrcorrection(result.Pvalue.to_list()))[1].astype(int)
        result["inter"]=inter_genes[topic]
        result["len_inter"]=len_inter_genes[topic]
        result["p_gene_topic"]=p_gene_topic[topic]
        result.sort_values(by="fdr",inplace=True, ascending=False)
        result.to_csv(f"{path_to_save}/{dataset} level {level} Enrichment Test {top_name}-topic {topic}.csv")
        
        
def topics_names(topics, enr_test_outcome, database):
    ''' 
    Procedure to assignt to each topic its most suitable name, that is the most suitable gene set from the database
    
    topics (list of strings) = list containing the names of the topic. They must be at least a subset of the "enr_test_outcome"'s keys
    enr_test_outcome (dictionary of dataframes) = each key correspond to a topic and each value to a a dataframe index by the gene sets (nameS)
    database (dictionary) = each key correspond to the name of the gene set and each value is a list of genes    
    '''
    
    topic_name={}
    for topic in topics:
        if len(enr_test_outcome[topic]) != 0:
            topic_name[topic]=[enr_test_outcome[topic].iloc[0].fdr.item(),enr_test_outcome[topic].index[0]]
        else:
            topic_name[topic]=["na","Not found"]
    topic_name=pd.DataFrame.from_dict(topic_name, orient="index")
    
    gain=len(set(topic_name[1]))
    
    while (len(set(topic_name[1])) != len(set(topic_name.index))):
        temp=topic_name.where(topic_name.duplicated(subset=1, keep=False)).dropna().sort_values(by=[1,0], ascending=False)
        temp["topic"]=temp.index
        temp=temp[temp.index!="Not found"]
        temp.set_index(1, inplace=True)
        for topic in list(set(temp[temp.index!="Not found"].index)):
            for top in temp.loc[topic][1:].topic:
                b=enr_test_outcome[top].drop(topic)
                for nuovo in b.index:
                    #Option 1: new name doesn't belong to any topic
                    if nuovo not in list(topic_name[1]):
                        topic_name[1].loc[top]=nuovo
                        topic_name[0].loc[top]=b.loc[nuovo].fdr
                        break
                    #Option 2: new name belongs to another topic --> check the -log10(fdr)
                    else:
                        f1=topic_name[topic_name[1]==nuovo][0].max()
                        f2=b.loc[nuovo].fdr
                        if f1>=f2: 
                            continue
                        else:
                            topic_name[1].loc[top]=nuovo
                            topic_name[0].loc[top]=b.loc[nuovo].fdr

        if len(set(topic_name[1])) > gain:

            gain=len(set(topic_name[1]))
        else:
            break
    print("Completeness:",len(set(topic_name[1])),"/",len(set(topic_name.index)))
    topic_name.columns=["fdr","name"]
    topic_name["inter"]=[enr_test_outcome[top].loc[topic_name.loc[top]["name"]].inter if topic_name.loc[top]["name"] != "Not found" else "Not found" for top in topic_name.index]
    topic_name["len_inter"]=[enr_test_outcome[top].loc[topic_name.loc[top]["name"]].len_inter if topic_name.loc[top]["name"] != "Not found" else "Not found" for top in topic_name.index]
    topic_name["len_gene_set"]=[len(database[l]) if l != "Not found" else "Not found" for l in topic_name["name"].to_list()]
    topic_name["p_gene_topic"]=[enr_test_outcome[top].loc[topic_name.loc[top]["name"]].p_gene_topic if topic_name.loc[top]["name"] != "Not found" else "Not found" for top in topic_name.index]

    return topic_name
   
   
   
####################################################################################################################################################################################################################### 
   
           

def col_clusters(lab, to_count, base):
    palette=[]
    for cla in sorted(set(lab[to_count])):
        counter=dict(collections.Counter(lab[lab[to_count]==cla][base]))
        inv_map = {v: k for k, v in counter.items()}
        palette.append(nmi.set_colors([f"{inv_map[list(inv_map.keys())[np.argmax(list(inv_map.keys()))]]}"]))
    return [item for sublist in palette for item in sublist]    
        

def hypergeom(my_genes,gene_set,M):
    n=len(gene_set)
    N=len(my_genes)
    x=len(list(set(gene_set) & set(my_genes)))
    return stats.hypergeom.sf(x-1, M, n, N)
   

def ens_to_name_l(l, trad):
    try:
        return list(trad.loc[l]["Gene name"])
    except:
        names=[]
        for ens in l:
            if ens in trad.index:
                names.append(trad.loc[ens]["Gene name"])
            else:
                names.append(ens)
        return names  


def intersection (lists):
    return list(set.intersection(*map(set,lists)))


def flat_list(t):
    return np.array(list(set([item for sublist in t for item in sublist])))


    
