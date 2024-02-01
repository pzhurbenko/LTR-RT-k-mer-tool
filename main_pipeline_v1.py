#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# маша задачи
1 Ускорить время подсчета, чтобы считалось не полтора часа
2 Базовый пайплайн с базовыми функциями


# In[1]:


import kmer_functions as kf
import time
import joblib


# In[4]:


all_seqs = "/home/pzhurb/Gypsy/Data/new_ann/Alo/data_good/LTR_seqs.fas"
rep_seqs = "/home/pzhurb/Gypsy/Data/new_ann/Alo/data_good/LTR_clusters.txt"
clades_annotated = "/home/pzhurb/Gypsy/Data/new_ann/Alo/Copia_clades.txt"
kmer_database = "/home/pzhurb/Gypsy/Data/new_ann/Alo/data_good/Copia_Gypsy/FTW_Copia_Gypsy.k13.fas"


# In[5]:


# Распарсить фасту на словарь
seqs_dict = kf.fasta_to_dict(all_seqs)
# Аннотированные названия клад и какие названия кластеров к ним относятся 
clades_dict = kf.parse_clades_file(clades_annotated)
# Распарсить файл с кластерами = Названия кластера: названия последовательностей в этом кластере
cluster_dict = kf.cluster_file_to_dict(rep_seqs)
# распарсить все кмеры на словарь (предварительно подготовить с помощью jellyfish)
kmer_database_dict = kf.fasta_to_dict_reverse(kmer_database)


# In[12]:


# Выбрать клады которые хочешь 
new_clades_dict = kf.extract_from_clade_file(clades_dict, \
    ['Angela'], clip = 3)
for i in new_clades_dict.keys():
    print(i)


# In[ ]:


# Если нужно объединить все выбранные клады, иначе пропускаем
joined_key = "|".join(list(set(new_clades_dict.keys())))
joined_value = []
joined_dict = {}
for i in new_clades_dict.values():
    joined_value.extend(i)
joined_dict[joined_key] = joined_value
new_clades_dict = joined_dict


# In[13]:


start_time = time.time()
k = 13
threashold = 30
fasta_out = 'clades_fasta_temp'
clades_uniq_nosingle_kmers = {}
stats_for_clades = {}
kmer_graph = {}
kmers_for_clades = {}
# Брать по одной кладе, и каждый кластер в кладе дополнить названиями последовательностей
yield_new_clades_dict = kf.yield_clades_dict(new_clades_dict, cluster_dict)
for temp_clade in yield_new_clades_dict:
    # temp_clade[0] - название клады
    # сколько кластеров в кладе
    len_temp_clade = len(temp_clade[1].values())
    print(temp_clade[0], len_temp_clade)
    # для текущей клады взять фасту в виде словаря
    clades_fasta = kf.clades_dict_to_clades_fasta(temp_clade[1], seqs_dict)
    # Записать фасту в файл, если это просто словарь
    kf.write_fasta_dict_to_fasta_file_2(clades_fasta, fasta_out + '.fas')
    # выполнить jellyfish
    kf.jellyfish_fasta(k, fasta_out)
    # Прочитать словарь к-меров после jellyfish
    kmers_dict = kf.fasta_to_dict_reverse(fasta_out + f".k{k}.fas")
    # отобрать уникальные к-меры
    uniq_kmers = kf.uniques_for_cluster(kmer_database_dict, kmers_dict)
    # Убрать сингтоны если кластер не синглтонный
    if len_temp_clade != 1:
        uniq_kmers = kf.filter_dict_by_value(uniq_kmers, '1')
        clades_uniq_nosingle_kmers[temp_clade[0]] = uniq_kmers
    else:
        clades_uniq_nosingle_kmers[temp_clade[0]] = uniq_kmers
        
    # граф в каких последовательностях и сколько раз встречается к-меры 
    temp_kmers_list = kf.find_uniq_kmers_in_cluster(clades_fasta, k, set(uniq_kmers), checkpoint = False)

    # Записать temp_kmers_list в словарь kmer_graph 
    kmer_graph[temp_clade[0]] = temp_kmers_list
    # Посчитать самые распространенные к-меры
    set_of_kmers, length, counts, dict_of_values = kf.abundant_kmers(temp_kmers_list, threashold)
    stats_for_clades[temp_clade[0]] = [length, counts, len(set_of_kmers), dict_of_values]
    # запишем к-меры
    kmers_for_clades[temp_clade[0]] = set_of_kmers
    
    end_time = time.time()
    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print("Elapsed time: ", round(elapsed_time, 1))

    # length = кол-во последовательностей, \
    # counts = минимальное кол-во последовательностей в которых встречается к-мер
    # для 'Angela,3,2,2,2,2' считалось 12 минут


# In[ ]:


# еще бы сделать фильтр для только уникальных против тех которые 2 раза встречаются


# In[15]:


kmers_for_clades


# In[23]:


joblib.dump(kmer_graph, 'kmer_graph_Copia_noAngela_4')


# In[23]:


def write_temp_kmers(temp_kmers_file, value):
    with open(temp_kmers_file, 'w') as file:
        for sub_value in value:
            file.write('>1' + '\n' + str(sub_value) + '\n')

def read_kmers_counts(temp_kmers_counts):
    kmers_counts_dict = {}
    with open(temp_kmers_counts, 'r') as file:
        for line in file:
            line = line.strip()  # Убираем символы новой строки и пробелы
            key, value = line.split()  # Разделяем строку по пробелу
            kmers_counts_dict[key] = int(value)
    return kmers_counts_dict


# In[24]:


kmers_for_clades_in_genomes = {}
list_of_jf = get_ipython().getoutput(' find /media/pzhurb/Diablo/Gypsy/genomes_kmers/Avena_jf/ -type f -name "*.jf"')
for key, value in kmers_for_clades.items():
    write_temp_kmers('temp_kmers.txt', value)
    for i in list_of_jf:
        get_ipython().system(' jellyfish query {i} -s temp_kmers.txt -o temp_kmers_counts.txt')
        total_kmers = get_ipython().getoutput(' jellyfish stats {i} | grep "Total" | awk \'{{print $$2}}\'')
        total_kmers = int(total_kmers[0])
        genome_name = i.split("/")[-1]
        kmers_counts = read_kmers_counts('temp_kmers_counts.txt') # словарь
        normalized_kmers_counts = {key: value / total_kmers for key, value in kmers_counts.items()}
        if key not in kmers_for_clades_in_genomes:
            kmers_for_clades_in_genomes[key] = [{genome_name:normalized_kmers_counts}]
        else:
            kmers_for_clades_in_genomes[key].extend([{genome_name:normalized_kmers_counts}])


# In[50]:


k = 13
threashold = 10
fas = '/home/pzhurb/Gypsy/Data/new_ann/Codes/kmer_prog_2.0/seqs_temp'
fas_dict = kf.fasta_to_dict(fas + '.fas')
kf.jellyfish_fasta(k, fas) # Исправить чтобы не учитывался .fas на конце а то получается .fas.fas
# Прочитать словарь к-меров после jellyfish
kmers_dict = kf.fasta_to_dict_reverse(fasta_out + f".k{k}.fas")


# In[54]:


temp_kmers_list = kf.find_uniq_kmers_in_cluster(fas_dict, k, set(kmers_dict), checkpoint = False)
set_of_kmers, length, counts, dict_of_values = kf.abundant_kmers(temp_kmers_list, threashold)


# In[82]:


seq_of_interest = 'chr3A:164931162-164936317'
for key, value in temp_kmers_list.items():
    if seq_of_interest in value[1]:
        print(value)


# In[ ]:




