#!/usr/bin/env python
# coding: utf-8

# In[1]:


from collections import Counter
import kcounter
from collections import defaultdict as dd
from Bio import SeqIO
import pandas as pd
import joblib
from copy import copy


# In[ ]:


def get_subdict_by_keys(dictionary, keys):
    subdict = {key: dictionary[key] for key in keys}
    return subdict


# In[ ]:


def get_subdict_by_values(dictionary, values):
    subdict = {k:v for k, v in dictionary.items() if v in values}
    return subdict


# In[ ]:


def common_in_two_dicts(dict1, dict2):
    common_pairs = {}
    for key in dict1:
        if key in dict2 and str(dict1[key]) == str(dict2[key]):
            common_pairs[key] = str(dict1[key])
    return common_pairs


# In[ ]:


def cluster_file_to_dict(file):
    df_clusters = pd.read_csv(file, sep = '\t', names = [0,1])
    cluster_dict = dd(list)
    for i in df_clusters.iterrows():
        cluster_dict[i[1][0]].append(i[1][1])
    return cluster_dict


# In[ ]:


def fasta_to_dict(fasta):
    seqs_dict = {}
    with open(fasta, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs_dict[record.id] = str(record.seq)
    return seqs_dict


# In[ ]:


def fasta_to_dict_reverse(fasta):
    seqs_dict = {}
    with open(fasta, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs_dict[str(record.seq)] = record.id
    return seqs_dict


# In[ ]:


def yield_seqs_by_cluster_name(cluster_dict, seqs_dict):
    for key, value in cluster_dict.items():
        the_cluster_dict = {key2:value2 for key2,value2 in seqs_dict.items() if key2 in value}
        yield [key, the_cluster_dict]


# In[ ]:


def filter_dict_by_value(input_dict, target_value):
    filtered_dict = {key: value for key, value in input_dict.items() if value != target_value}
    return filtered_dict


# In[ ]:


# def find_uniq_kmers(seq, kmer, count_kmer_only_once = True):
#     kmer_dict = kcounter.count_kmers(seq, kmer, canonical_kmers=True)
#     if count_kmer_only_once == True:
#         kmer_dict.update((k,1) for k in kmer_dict)
#         return kmer_dict
#     else:
#         kmer_dict.update((k, int(v)) for k, v in kmer_dict.items()) # избавимся от float
#         return kmer_dict


# In[ ]:


# def find_uniq_kmers_in_cluster(the_cluster_dict, kmer, count_kmer_only_once = True, remove_singletones = False):
#     kmer_dict = {}
#     for key, value in the_cluster_dict.items():
#         kmer_dict = dict(Counter(kmer_dict) +\
#             Counter(find_uniq_kmers(value, kmer, count_kmer_only_once = count_kmer_only_once)))
#     if remove_singletones != False:
#         kmer_dict = filter_dict_by_value(kmer_dict, 1)
#     return kmer_dict


# In[ ]:


# def filter_kmers(all_clusters_and_kmers, uniq_kmers):
#     for i in range(0, len(all_clusters_and_kmers)):
#         for j in list(all_clusters_and_kmers[i][1]): # list чтобы не было ошибки при итерировании
#             if j not in uniq_kmers:
#                 del all_clusters_and_kmers[i][1][j]
#     return all_clusters_and_kmers # измененный. Чтобы сохранить оригинал, надо сделать копию


# In[ ]:


# def count_kmers_in_clusters(the_cluster_dict, how_many_clusters_to_parse, count_kmer_only_once = True):
#     all_clusters_and_kmers = []
#     all_cluster_size = []
#     for _num_ in range(0, how_many_clusters_to_parse):
#         try:
#             i = next(the_cluster_dict)
#             all_cluster_size.append(len(i[1]))
#             if len(i[1]) == 1:
#                 for key, value in i[1].items():
#                     all_clusters_and_kmers.append([i[0], find_uniq_kmers(value, 13, count_kmer_only_once = count_kmer_only_once)])
#             else:
#                 all_clusters_and_kmers.append([i[0], find_uniq_kmers_in_cluster(i[1], 13, count_kmer_only_once = count_kmer_only_once)])
#         except StopIteration:
#             return all_clusters_and_kmers
#     return all_clusters_and_kmers


# In[ ]:


# Если в базе всех к-меров и к-мерах кластера совпдают числа, то это уникальные для кластера

def uniques_for_cluster(kmer_database_dict, cluster_kmers_dict):
    common_kmer_seqs = get_subdict_by_keys(kmer_database_dict, list(cluster_kmers_dict))
    uniq_kmers_for_cluster = common_in_two_dicts(common_kmer_seqs, cluster_kmers_dict)
    return uniq_kmers_for_cluster


# In[ ]:


def write_fasta_dict_to_fasta_file(dictionary, fasta_output):
    with open (fasta_output, "w") as f:
        for i in dictionary:
            for key, value in i[1].items():
#                f.write(">" + key + "|" + i[0] + "\n" + value + "\n") # Записать с номерами кластеров через "|"
                f.write(">" + key + "\n" + value + "\n")


# In[ ]:


# из файла с кладами (Copia_clades.txt) формата chr1A:100707300-100718043\tCopia\tSIRE\t1 создать словарь 
# где будут названиями клад:[названия последовательностей]. Берутся самые мелкие клады

def parse_clades_file(file):
    records = []
    with open(file, 'r') as handle:
        for line in handle:
            line = line.rstrip()
            records.append(line.split("\t"))
    set_of_clades = set()
    for clust in records:
        set_of_clades.update([(',').join(clust[1:])])
    dict_of_clades = {}
    for clade in set_of_clades:
        for clust in records:
            cluster_name = (',').join(clust[1:])
            if clade == cluster_name:
                if clade in dict_of_clades:
                    dict_of_clades[cluster_name].append(clust[0])
                else:
                    dict_of_clades[cluster_name] = [clust[0]]
    return dict_of_clades


# In[ ]:


def yield_clades_dict(clades_dict, cluster_dict):
    for key, value in clades_dict.items():
        temp_clade = [key, get_subdict_by_keys(cluster_dict, value)]
        yield temp_clade


# In[ ]:


def jellyfish_fasta(k, fasta_out):
    get_ipython().system(' jellyfish count -m {k} -s 1G -t 7 -C -o {fasta_out}.k13.jr {fasta_out}.fas')
    get_ipython().system(' jellyfish dump {fasta_out}.k13.jr > {fasta_out}.k13.fas')


# In[ ]:


# Записать фасту в файл, если это просто словарь, а не generator object 
def write_fasta_dict_to_fasta_file_2(dictionary, fasta_output):
    with open (fasta_output, "w") as f:
        for key, value in dictionary.items():
            f.write(">" + key + "\n" + value + "\n")
    #       f.write(">" + key + "|" + i[0] + "\n" + value + "\n") # Записать с номерами кластеров через "|"


# In[ ]:


# Получить один фаста - словарь из словаря new_clades_dict (а не парсить по одной кладе)
def clades_dict_to_clades_fasta(clades_dict, seqs_dict):
    seqs_names = []
    for i in clades_dict.values():
        for j in i:
            seqs_names.append(j)
    clades_fasta = get_subdict_by_keys(seqs_dict, seqs_names)
    return clades_fasta


# In[ ]:


def get_seqs_by_cluster_name(cluster_dict, seqs_dict):
    for key, value in cluster_dict.items():
        the_cluster_dict = {key2:value2 for key2,value2 in seqs_dict.items() if key2 in value}
        yield [key, the_cluster_dict]


# In[ ]:


def extract_from_clade_file(clades_dict, feature_list, clip = 0):
    new_clades_dict = {}
    if clip == 0:
        for key, value in clades_dict.items():
            for feature in feature_list:
                if feature in key:
                    new_clades_dict[key] = value
    else:
        for key, value in clades_dict.items():
            for feature in feature_list:
                if feature in key:
                    key = key.split(',')
                    key = key[:clip]
                    key = ','.join(key)
                    if key not in new_clades_dict:
                        new_clades_dict[key] = copy(value) # ваще не понимаю про словари
                    else:
                        new_clades_dict[key].extend(value)
    return new_clades_dict


# In[ ]:





# In[ ]:





# In[ ]:


def merge_dictionaries(dict1, dict2):
    result = {}
    for key in dict1:
        if key in dict2:
            result[key] = [dict1[key], dict2[key]]
    return result

# {b': 2, 'c': 3} + {'a': 1, b': 4, 'c': 5} = {'b': [2, 4], 'c': [3, 5]}


# In[ ]:


def merge_merged_dictionaries(dict1, dict2):
    result = {}
    union_set = set().union(dict1, dict2)
    intersection_set = set(dict1).intersection(dict2)
    for key in union_set:
        if key in intersection_set:
            merged_value = [0] * len(dict1[key][0])
            for i in range(0, len(dict1[key][0])):
                merged_value[i] = dict1[key][0][i] + dict2[key][0][i]
            merged_value = [merged_value, dict1[key][1] + dict2[key][1]]
            result[key] = merged_value
    for key in set(dict1).difference(dict2):
        result[key] = dict1[key]
    for key in set(dict2).difference(dict1):
        result[key] = dict2[key]    
    return result

# merge_merged_dictionaries({'kmer': [[1, 1], [clust_name1]]}, {'kmer': [[1, 1], [clust_name2]]}) =\
#    {'kmer': [[2, 2], [clust_name1, clust_name2]]}


# In[ ]:


def filter_dict(input_dict, filter_keys):
    filtered_dict = {key: input_dict[key] for key in filter_keys if key in input_dict}
    return filtered_dict
# {'a': 1, 'b': 2, 'c': 3, 'd': 4} + {'a', 'c', 'e'} = {'c': 3, 'a': 1}


# In[ ]:


def find_uniq_kmers(seq, k, seq_name = None, good_kmers = None):
    kmer_dict = kcounter.count_kmers(seq, k, canonical_kmers=True)
    if good_kmers != None:
        filtered_dict = filter_dict(kmer_dict, good_kmers)
    else:
        filtered_dict = kmer_dict
    filtered_dict.update((key, int(v)) for key, v in filtered_dict.items()) # избавимся от float
    kmer_dict_uniq = filtered_dict.copy()
    kmer_dict_uniq.update((key, 1) for key in filtered_dict)
    merged_dict = merge_dictionaries(filtered_dict, kmer_dict_uniq)
    if seq_name != None:
        for key, value in merged_dict.items():
            new_value = [value, [seq_name]]
            merged_dict[key] = new_value
        return merged_dict
    else:
        return merged_dict

# {'AGGACTAGAGTTG': [1, 2] that stands for [uniq count, real count]


# In[ ]:


def find_uniq_kmers_in_cluster(seq_dict, k, good_kmers = None, checkpoint = False):
    kmer_dict = {}
    count = 0
    for key, value in seq_dict.items():
        if checkpoint != False:
            count += 1
            if count // 50 == 0:
                print(len(seq_dict), count)
        temp_dict = find_uniq_kmers(value, k, seq_name = key, good_kmers = good_kmers)
        kmer_dict = merge_merged_dictionaries(kmer_dict, temp_dict)
    return kmer_dict

# seq_dict = {'seq_name': 'seq', 'seq_name2': 'seq2', ....}


# In[2]:


def abundant_kmers(kmers_counts_names, threshold):
    data_dict = kmers_counts_names.copy()
    data_list = [[key] + value for key, value in data_dict.items()]
    sorted_data = sorted(data_list, key=lambda x: x[1][1], reverse=True)
    set_all_seq_names = set()
    for i in data_list:
        set_all_seq_names.update(i[2])
    set_of_kmers = set()
    dict_of_clust_names = {}
    for i in sorted_data:
        # смысл избежать избыточного добавления к-меров если они уже встречаются во всех последовательностях, \
        # для которые уже были найдены к-меры в числе бОльшем чем treashold
        # Проверяем, что все clust_names из i присутствуют в словаре dict_of_clust_names
        all_keys_present = all(key in dict_of_clust_names for key in i[2])
        if all_keys_present:
            cluster_counts = [dict_of_clust_names[key] for key in i[2]]
            if min(cluster_counts) >= threshold:
                continue
        set_of_kmers.update([i[0]])
        for j in i[2]: # clust_names
            if j not in dict_of_clust_names:
                dict_of_clust_names[j] = 1
            else:
                dict_of_clust_names[j] += 1

        if set(dict_of_clust_names) == set_all_seq_names: #Узнаем минимальное значение в словаре, когда он заполнится
            min_value = min(dict_of_clust_names, key=dict_of_clust_names.get)
            if dict_of_clust_names[min_value] >= threshold:
                break
    return set_of_kmers, len(set(dict_of_clust_names)), i[1][1], dict_of_clust_names
    
# kmers_counts_names = {'kmer1': [[4, 2], [clust_name1, clust_name2]], 'kmer2': [[2, 2], [clust_name1, clust_name2]]}
# минимальное количество кмеров = threshold в каждом кластере. 
# Выдвется сколько сиквенсов в кластере (length = len(set(dict_of_clust_names)))
# и в скольких кластерах стречается самый малопредставленный кмер  (counts = i[1][1])


# In[ ]:





# In[ ]:




