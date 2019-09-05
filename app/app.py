from flask import Flask, render_template
import os
import pandas as pd
from pandas.compat import StringIO
import numpy as np
from numpy import loadtxt
import sys
import json
from pprint import pprint
import objectpath
import csv
import re
import matplotlib.pyplot as plt
import json
import io
import http.client
import requests
from pprint import pprint
import itertools
import scipy as sp
from scipy.spatial import distance
from sklearn.metrics.pairwise import pairwise_distances
from clustergrammer import Network
#from clustergrammer_widget import *
from flask import Flask, render_template, request, redirect, jsonify, send_from_directory, abort, Response, send_file
import urllib.request
import h5py
import h5sparse
from scipy.sparse import csc_matrix, csr_matrix, coo_matrix, vstack
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
# os.chdir('./scripts')
from .gvm_scripts import *
from .vae_scripts import *
import keras
from keras import backend as K
import tensorflow

app = Flask(__name__, static_folder="./static/dist", template_folder="./static")

@app.route("/")
def index():
    DrOI = request.args.get('drug')
    print(DrOI)
    # DrOI = "acarbose"
    DrugMatrix ={}
    for line in urllib.request.urlopen('http://amp.pharm.mssm.edu/lincs-playground/index.php/s/bl6rC50ALsbXfw1/download'):
        label, genelist = line.decode().split('\t\t', maxsplit=1) 
        genelist_split = genelist.strip().split('\t')
        DrugMatrix[label] = genelist_split
    DrugMatrix = {x.replace('/', '_'): v  
        for x, v in DrugMatrix.items()}
    ### DrugMatrix drug input
    ## generate a list of searchable keys to reduce dict
    DrugMatrix_keys = pd.DataFrame(list(DrugMatrix.keys()))
    DrugMatrix_keys.columns = ["sigs"]
    Drug_Matrix_DrOI = DrugMatrix_keys[DrugMatrix_keys["sigs"].apply(lambda s: bool(re.compile(str(DrOI), re.IGNORECASE).search(str(s))))]
    ## reduce dict
    Drug_matrix_sigs_reduced = list(Drug_Matrix_DrOI["sigs"])
    #DrugMatrix_sigs = {k: DrugMatrix[k] for k in list(Drug_matrix_sigs_reduced["sigs"])} # total sigs
    ## up sigs
    Drug_matrix_up_sigs_reduced = Drug_Matrix_DrOI[Drug_Matrix_DrOI["sigs"].apply(lambda s: bool(re.compile(str("-up"), re.IGNORECASE).search(str(s))))]
    DrugMatrix_up_sigs= {k: DrugMatrix[k] for k in list(Drug_matrix_up_sigs_reduced["sigs"])}
    for a in list(Drug_matrix_up_sigs_reduced["sigs"]):
        DrugMatrix_up_sigs_save = DrugMatrix_up_sigs[a]
        print(a)
        #with open(a + "_DrugMatrix_up_sig.json", "w") as f:
            #json.dump(DrugMatrix_up_sigs_save, f)
            
    ## down sigs
    Drug_matrix_down_sigs_reduced = Drug_Matrix_DrOI[Drug_Matrix_DrOI["sigs"].apply(lambda s: bool(re.compile(str("-dn"), re.IGNORECASE).search(str(s))))]
    DrugMatrix_down_sigs= {k: DrugMatrix[k] for k in list(Drug_matrix_down_sigs_reduced["sigs"])}
    for b in list(Drug_matrix_down_sigs_reduced["sigs"]):
        DrugMatrix_down_sigs_save = DrugMatrix_down_sigs[b]
        print(b)
        #with open(b + "_DrugMatrix_down_sig.json", "w") as f:
            #json.dump(DrugMatrix_down_sigs_save, f)
            
            
    ### L1000 DRUG CARD DRUG INPUT
    DrOI_df = metadata[metadata["pert_desc"] == DrOI]
    DrOI_pert_ids = list(DrOI_df["pert_id"])
    DrOI_up_signatures = {k: L1000_up_lookup[k] for k in (DrOI_pert_ids)}
    DrOI_up_no_perts = {k: v for d in DrOI_up_signatures.values() for k, v in d.items()}
    DrOI_up_drug_sigs = list(DrOI_up_no_perts.keys())
    DrOI_down_signatures = {k: L1000_down_lookup[k] for k in (DrOI_pert_ids)}
    DrOI_down_no_perts = {k: v for d in DrOI_down_signatures.values() for k, v in d.items()}
    DrOI_down_drug_sigs = list(DrOI_down_no_perts.keys())
    DrOI_all_sigs = set(DrOI_up_drug_sigs) & set (DrOI_down_drug_sigs)
    DrOI_all_sigs_up = [s + "_up" for s in DrOI_all_sigs]
    DrOI_all_sigs_down = [s + "_down" for s in DrOI_all_sigs]
    ######## NEW CODE
    DrOI_all_sigs_display = [DrOI + "_" + s for s in list(DrOI_all_sigs)]
    ########
    for a in DrOI_all_sigs:
        L1000_up_json_file = DrOI_up_no_perts[a]
        L1000_down_json_file = DrOI_down_no_perts[a]
        print(a)
        #with open(a + "_L1000_up_sig.json", "w") as f:
            #json.dump(L1000_up_json_file, f)
        #with open(a + "_L1000_down_sig.json", "w") as f:
            #json.dump(L1000_down_json_file, f)
            
    ### CREEDS DRUG CARD 
    #for a in loop_iteration:
    CREEDS_URL = 'http://amp.pharm.mssm.edu/CREEDS/'
    CREEEDS_Drug_response = requests.get(CREEDS_URL + 'search', params={'q':DrOI})
    if CREEEDS_Drug_response.status_code == 200:
        #pprint(CREEEDS_Drug_response.json())
        #json.dump(CREEEDS_Drug_response.json(), open(DrOI + '_api1_result.json', 'w'), indent=4)
        CREEDS_drug_output_df = pd.DataFrame(CREEEDS_Drug_response.json())
        CREEDS_drug_output_ids = list(CREEDS_drug_output_df["id"])
        
        CREEDS_drug_output_ids_up = ["CREEDS_" + s + "_up" for s in CREEDS_drug_output_ids]
        CREEDS_drug_output_ids_down = ["CREEDS_" + s + "_down" for s in CREEDS_drug_output_ids]
        
        CREEDS_all_down_genes = []
        CREEDS_all_up_genes = []
        CREEDS_desc = []
        for a in CREEDS_drug_output_ids:
            CREEDS_URL = 'http://amp.pharm.mssm.edu/CREEDS/'
            CREEDS_drug_sigs_response = requests.get(CREEDS_URL + 'api', params={'id':str(a)})
            if CREEDS_drug_sigs_response.status_code == 200:
                CREEDS_drug_sigs_response_json = CREEDS_drug_sigs_response.json()
                
                ## up genes
                CREEDS_drug_sigs_up_genes = CREEDS_drug_sigs_response_json['up_genes']
                CREEDS_drug_sigs_up_genes_df = pd.DataFrame(CREEDS_drug_sigs_up_genes) # this is the up genes dataframe
                CREEDS_drug_sigs_up_genes_df.columns = ["Genes", "Score"]
                filename1 = (a + "_CREEDS_drug_sig_up_genes.csv")
                #CREEDS_drug_sigs_up_genes_df.to_csv(filename1) # this saves the df as a csv
                desc = (a + "_" + DrOI + "_" + CREEDS_drug_sigs_response_json["geo_id"])
                CREEDS_desc.append(desc)
                CREEDS_all_up_genes.append(list(CREEDS_drug_sigs_up_genes_df["Genes"]))
                
                ## down genes
                CREEDS_drug_sigs_down_genes = CREEDS_drug_sigs_response_json['down_genes']
                CREEDS_drug_sigs_down_genes_df = pd.DataFrame(CREEDS_drug_sigs_down_genes)# this is the down genes dataframe
                CREEDS_drug_sigs_down_genes_df.columns = ["Genes", "Score"]
                filename2 = (a + "_CREEDS_drug_sig_down_genes.csv")
                CREEDS_all_down_genes.append(list(CREEDS_drug_sigs_down_genes_df["Genes"]))
                #CREEDS_drug_sigs_down_genes_df.to_csv(filename2)
                #CREEDS_drug_sigs_down_genes = CREEDS_drug_sigs_response_json['down_genes'] # this saves the df as a csv
                print(filename2)
                
                up_keys = ['up_genes']
                gene_dict_up = {x:CREEDS_drug_sigs_response_json[x] for x in up_keys}
                gene_dict_up = {"CREEDS_" + a + "_" + k: v for k, v in gene_dict_up.items()}
                
                down_keys = ['down_genes']
                gene_dict_down = {x:CREEDS_drug_sigs_response_json[x] for x in down_keys}
                gene_dict_down = {"CREEDS_" + a + "_" + k: v for k, v in gene_dict_down.items()}
                
#storing genes in gene sets and signatures 
    ### CREEDS DISEASE CARD (DRUG INPUT)
    # RETURNS THE do_id, geo_id, and disease name in a dictionary
    CREEDS_GSE = {
        row['id']: [row['geo_id'], row["disease_name"]]
        for row in CREEDS_data
    }
    ## filter by DrOI need icd9 codes for proper conversion and query through CREEDS
    droi_search =EMR_data_df[EMR_data_df['Drug_Name'].apply(lambda s: bool(re.compile(DrOI, re.IGNORECASE).search(s)))]
    droi_search_top5 = droi_search[0:10]
    EMR_top_disease_from_drug = droi_search_top5["ICD9"]
    #top_disease_from_drug = EMR_top_disease_from_drug[0:5]
    ## build a datatable of all the ICD-9 CM diagnosis codes families (i.e no decimal points)
    EMR_top_disease_from_drug_df = pd.DataFrame(EMR_top_disease_from_drug, columns=['ICD9'])
    EMR_top_disease_from_drug_df['ICD9_wildcard'] = EMR_top_disease_from_drug_df['ICD9'].apply(lambda code: code.split('.')[0])
    #EMR_top_disease_from_drug_df.head()
    icd9_to_doid_final['ICD9_wildcard'] = icd9_to_doid_final['ICD9'].apply(lambda code: str(code).split('.')[0])
    #icd9_to_doid_final.head()
    df_joined = pd.merge(
        left=EMR_top_disease_from_drug_df, left_on='ICD9_wildcard',
        right=icd9_to_doid_final, right_on='ICD9_wildcard',
        how='inner',
        suffixes=(
            '_left',
            '_right',
        )
    )
    CREEDS_drug_ids = pd.DataFrame(set(df_joined.CREEDS_drug_id))
    CREEDS_drug_ids_list = list(set(df_joined.CREEDS_drug_id))
    #CREEDS_GSE.keys()
    #CREEDS_drug_ids_list
    CREEDS_Drug_Final = dict((k, CREEDS_GSE[k]) for k in CREEDS_drug_ids_list)
    CREEDS_drug_final_df = pd.DataFrame(CREEDS_Drug_Final).T
    CREEDS_drug_final_df.columns = ["GSE_ID", "DISEASE"]
    #CREEDS_drug_final_df # DISPLAY THIS DATAFRAME
    ### CREEDS DISEASE CARD FROM DRUG INPUT API
    CREEDS_drug_final_diseases = CREEDS_drug_final_df.DISEASE
    CREEDS_drug_final_GSE_ID = CREEDS_drug_final_df.GSE_ID
    ## CREEDS DISEASE CARD FROM DISEASE QUERY 
    CREEDS_disease_output_ids_up = ["CREEDS_" + s + "_up" for s in CREEDS_drug_final_diseases]
    CREEDS_disease_output_ids_down = ["CREEDS_" + s + "_down" for s in CREEDS_drug_final_diseases]
    loop_iteration = np.arange(0, len(CREEDS_drug_final_diseases))
    loop_iteration = list(loop_iteration)
    CREEDS_total_api_df = []
    CREEDS_all_disease_up_genes = []
    CREEDS_all_disease_down_genes = []
    for a in loop_iteration:
        CREEDS_URL = 'http://amp.pharm.mssm.edu/CREEDS/'
        CREEEDS_Disease_response = requests.get(CREEDS_URL + 'api', params={'id':CREEDS_drug_ids_list[a]})
        if CREEEDS_Disease_response.status_code == 200:
            #pprint(CREEEDS_Disease_response.json())
            #json.dump(CREEEDS_Drug_response.json(), open(CREEDS_drug_final_GSE_ID[a] + '_api1_result.json', 'w'), indent=4)
            CREEEDS_Disease_response_json = CREEEDS_Disease_response.json()
            
            ## up genes
            CREEDS_disease_sigs_up_genes = CREEEDS_Disease_response_json['up_genes']
            CREEDS_disease_sigs_up_genes_df = pd.DataFrame(CREEDS_disease_sigs_up_genes) # this is the up genes dataframe
            CREEDS_disease_sigs_up_genes_df.columns = ["Genes", "Score"]
            #desc = (a + "_" + DrOI + "_" + CREEEDS_Disease_response["geo_id"])
            #CREEDS_desc.append(desc)
            CREEDS_all_disease_up_genes.append(list(CREEDS_disease_sigs_up_genes_df["Genes"]))
            
            filename1 = (str(CREEDS_drug_ids_list[a]) + "_CREEDS_disease_sig_up_genes.csv")
            #CREEDS_disease_sigs_up_genes_df.to_csv(filename1) # this saves the df as a csv
            
            
            ## down genes
            CREEDS_disease_sigs_down_genes = CREEEDS_Disease_response_json['down_genes']
            CREEDS_disease_sigs_down_genes_df = pd.DataFrame(CREEDS_disease_sigs_down_genes) # this is the up genes dataframe
            CREEDS_disease_sigs_down_genes_df.columns = ["Genes", "Score"]
            CREEDS_all_disease_down_genes.append(list(CREEDS_disease_sigs_down_genes_df["Genes"]))
            
            filename2 = (str(CREEDS_drug_ids_list[a]) + "_CREEDS_disease_sig_down_genes.csv")
            #CREEDS_disease_sigs_down_genes_df.to_csv(filename2) # this saves the df as a csv
            print(filename2)
            # entire json
            #json.dump(response.json(), open(a + '_CREEDS_Disease_sig.json', 'w'), indent=4) # if the user wants the entire json, they can download this
                
            
            #CREEEDS_Drug_response_df = pd.DataFrame(CREEEDS_Drug_response_json)
            #CREEEDS_Drug_response_df # This will be the dataframe to return
            #CREEDS_total_api_df.append(CREEEDS_Drug_response_df)
    #CREEDS_total_api_df = pd.concat(CREEDS_total_api_df, axis =1)
    #CREEDS_total_api_df.T ## display this datatable
    ### GENESHOT API and further integration
    GENESHOT_URL = 'http://amp.pharm.mssm.edu/geneshot/api'
    query_string = '/search/%s'
    search_term = DrOI
    # true query from geneshot
    response = requests.get(
        GENESHOT_URL + query_string % (search_term)
    )
    if not response.ok:
        raise Exception('Error during query')
    data = json.loads(response.text)
    #print(data)
    ## GENESHOT QUERY USING AutoRIF
    GENESHOT_URL = 'http://amp.pharm.mssm.edu/geneshot/api'
    query_string = '/search/auto/%s'
    search_term = 'wound healing' # this will be the user input 
    geneshot_response = requests.get(
        GENESHOT_URL + query_string % (search_term)
    )
    if not geneshot_response.ok:
        raise Exception('Error during query')
    geneshot_data = json.loads(geneshot_response.text)
    #print(geneshot_data)
    geneshot_gene_df = geneshot_data["gene_count"]
    geneshot_gene_list = list(geneshot_gene_df.keys()) # this extracts the genes from the json. We can then resend this through the geneshot api
    geneshot_gene_list_commas = ",".join(geneshot_gene_list) # can save this as a csv. 
    geneshot_gene_df1 = pd.DataFrame(geneshot_gene_df).T
    geneshot_gene_df1.columns = ["Pubmed Count", "Publication Count/Total Publications"]
    #write the geneshot pubmed data
    #geneshot_gene_df1.to_csv(search_term + "_geneshot_pubmed_counts.csv")
    query_string = '/associate/%s/%s'
    similarity_matrix = 'coexpression' # we can make this dynamic. Parameters: (generif, tagger, autorif, coexpression, enrichr)
    gene_symbols = geneshot_gene_list_commas
    coexpression_response = requests.get(
        GENESHOT_URL + query_string % (similarity_matrix, gene_symbols)
    )
    if not coexpression_response.ok:
        raise Exception('Error during query')
    coexpression_data = json.loads(coexpression_response.text) # this will be the coexpression json they can download
    geneshot_coexp_ass = {"GENESHOT_coexpression":list(coexpression_data["association"].keys())}
    query_string = '/associate/%s/%s'
    similarity_matrix = 'generif' # we can make this dynamic. Parameters: (generif, tagger, autorif, coexpression, enrichr)
    gene_symbols = geneshot_gene_list_commas
    generif_response = requests.get(
        GENESHOT_URL + query_string % (similarity_matrix, gene_symbols)
    )
    if not generif_response.ok:
        raise Exception('Error during query')
    generif_data = json.loads(generif_response.text) # this will be the coexpression json they can download
    geneshot_generif = {"GENESHOT_generif":list(generif_data["association"].keys())}
    query_string = '/associate/%s/%s'
    similarity_matrix = 'tagger' # we can make this dynamic. Parameters: (generif, tagger, autorif, coexpression, enrichr)
    gene_symbols = geneshot_gene_list_commas
    tagger_response = requests.get(
        GENESHOT_URL + query_string % (similarity_matrix, gene_symbols)
    )
    if not tagger_response.ok:
        raise Exception('Error during query')
    tagger_data = json.loads(tagger_response.text) # this will be the coexpression json they can download
    geneshot_tagger = {"GENESHOT_tagger":list(tagger_data["association"].keys())}
    query_string = '/associate/%s/%s'
    similarity_matrix = 'tagger' # we can make this dynamic. Parameters: (generif, tagger, autorif, coexpression, enrichr)
    gene_symbols = geneshot_gene_list_commas
    autorif_response = requests.get(
        GENESHOT_URL + query_string % (similarity_matrix, gene_symbols)
    )
    if not autorif_response.ok:
        raise Exception('Error during query')
    autorif_data = json.loads(autorif_response.text) # this will be the coexpression json they can download
    geneshot_autorif = {"GENESHOT_autorif":list(autorif_data["association"].keys())}
    query_string = '/associate/%s/%s'
    similarity_matrix = 'tagger' # we can make this dynamic. Parameters: (generif, tagger, autorif, coexpression, enrichr)
    gene_symbols = geneshot_gene_list_commas
    enrichr_response = requests.get(
        GENESHOT_URL + query_string % (similarity_matrix, gene_symbols)
    )
    if not enrichr_response.ok:
        raise Exception('Error during query')
    enrichr_data = json.loads(enrichr_response.text) # this will be the coexpression json they can download
    geneshot_enrichr = {"GENESHOT_enrichr":list(enrichr_data["association"].keys())}
    #### GMT formation from these datasets (NO X2K)
    #### format = TITLE \t\ Description \t\ Genes
    def merge_two_dicts(x, y):
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z
    DrOI_up_no_perts_cp = {"L1000FWD_" + k +"_up": v for k, v in DrOI_up_no_perts.items()}
    DrOI_down_no_perts_cp = {"L1000FWD_" + k +"_down": v for k, v in DrOI_down_no_perts.items()}
    ## Genes
    DrugMatrix_Drug_Genes = merge_two_dicts(DrugMatrix_up_sigs , DrugMatrix_down_sigs)
    DrugMatrix_Drug_Genes = {"DrugMatrix_" + k: v for k, v in DrugMatrix_Drug_Genes.items()}
    L1000_Drug_Genes = merge_two_dicts(DrOI_up_no_perts_cp,DrOI_down_no_perts_cp)
    CREEDS_up_Genes = {
        CREEDS_drug_output_ids_up[a]: CREEDS_all_up_genes[a]
        for a in range(len(CREEDS_drug_output_ids_up))
    }
    CREEDS_down_Genes = {
        CREEDS_drug_output_ids_down[a]: CREEDS_all_down_genes[a]
        for a in range(len(CREEDS_drug_output_ids_up))
    }
    CREEDS_Disease_up_Genes = {
        CREEDS_disease_output_ids_up[a]: CREEDS_all_disease_up_genes[a]
        for a in range(len(CREEDS_drug_output_ids_up))
    }
    CREEDS_Disease_down_Genes = {
        CREEDS_disease_output_ids_down[a]: CREEDS_all_disease_down_genes[a]
        for a in range(len(CREEDS_drug_output_ids_up))
    }
    total_genes = merge_two_dicts(L1000_Drug_Genes, DrugMatrix_Drug_Genes)
    total_genes =merge_two_dicts(total_genes, geneshot_coexp_ass)
    total_genes =merge_two_dicts(total_genes, geneshot_generif)
    total_genes =merge_two_dicts(total_genes, geneshot_tagger)
    total_genes =merge_two_dicts(total_genes, geneshot_autorif)
    total_genes =merge_two_dicts(total_genes, geneshot_enrichr)
    total_genes =merge_two_dicts(total_genes, CREEDS_up_Genes)
    total_genes =merge_two_dicts(total_genes, CREEDS_down_Genes)
    total_genes =merge_two_dicts(total_genes, CREEDS_Disease_up_Genes)
    total_genes =merge_two_dicts(total_genes, CREEDS_Disease_down_Genes)
    ## you will need to change this path file. 
    print("./data/" + DrOI)
    with open ("./data/" + DrOI+".gmt", "w") as file:
        for k in list(total_genes.keys()):
            file.write(k + '\t')
            #file.write('\t'+'na')
            file.write("\t".join(total_genes[k]))
            file.write('\n')

            

#AUTO ENCODER CODE

    

    all_genes = pd.read_csv('./data/AE/genes_info.tsv', sep='\t', index_col=0)
    gmt_fname = "./data/" + DrOI + '.gmt'
    lib_name = os.path.splitext(gmt_fname.rsplit('/', 1)[-1])[0]
    gvm_fname = './data/' + lib_name + '.h5'
    print(gvm_fname)
    formatted_gvm_fname = './data/' + lib_name + '_FORMATTED.h5'
    
    if os.path.isfile(gvm_fname): 
        gvm = open_gvm(gvm_fname)
    else:
        gvm = convert_genesetlist(get_genesetlist(gmt_fname, 'gmt_fname'), to='gvm_h5', output_fname=gvm_fname)
    
    summary = format_gvm_h5(gvm_fname = gvm_fname, all_genes = all_genes,
                       output_fname = formatted_gvm_fname, max_gs_loss=1.0, min_gs_size=1,
                       overwrite = False, return_value='summary')

    n_labels, n_genes = get_gvm_size(formatted_gvm_fname)
    (n_labels, n_genes)

    group = 'AE' # vanilla autoencoder

    batch_size = 128
    m = 1000 # middle dimension
    l = 50 # latent dimension

    model = build_vae(input_dim=n_genes, middle_dim = m, latent_dim = l, 
                    batch_size=batch_size, optimizer='Adamax', lr=.001)
    vae, enc, dec = (model['vae'], model['enc'], model['dec'])
    vae.load_weights('./models/%s/weights/%04dm_%04dl.h5'%(group, m, l))     

    z = enc.predict_generator(
    GeneVec_Generator(formatted_gvm_fname, gvm_path='gvm', batch_size=1000, shuffle=False),
    workers=4, use_multiprocessing=True, verbose=0)

    euc_dist = pairwise_distances(z, metric='euclidean')
    cos_sim = cosine_similarity(z)
    labels = open_gvm(formatted_gvm_fname)['idx']

    euc_dist_df = pd.DataFrame(euc_dist, index=labels, columns=labels)
    cos_sim_df = pd.DataFrame(cos_sim, index=labels, columns=labels) 
    euc_dist_df.iloc[:5, :5]
    cos_sim_df.iloc[:5, :5]
    euc_dist_df.to_pickle('./data/%s_DIST_EUC.pkl'%lib_name)
    cos_sim_df.to_pickle('./data/%s_DIST_COS.pkl'%lib_name)
    cos_sim_df2 = pd.read_pickle('./data/%s_DIST_COS.pkl'%lib_name)

    return render_template("index.html")


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=3005, debug=True)