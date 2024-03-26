import glob
import pandas as pd
import os
import numpy as np
from gensim.models.word2vec import Word2Vec
from gensim.models.keyedvectors import KeyedVectors
import gensim
import multiprocessing


def get_aa_mutations(tsv_path):        # Creating the dictionary

    aa_mutations = []
    file = tsv_path
    df = pd.DataFrame()
    for file in glob.glob(files):      # Merging tsv files downloaded from GISAID
        temp_df = pd.read_csv(file,sep='\t')
        df = df.append(temp_df)
    df = df.drop_duplicates(subset='seqName')

    df['epi_id'] = pd.Series(dtype='str')
    df.drop('index', axis=1, inplace=True)
    df = df.reset_index()

    for i,row in df.iterrows():
        temp_name = row['seqName']
        temp_name = temp_name.split('|')
        temp_name = temp_name[1]
        df.at[i,'epi_id'] = temp_name

    df[['aaSubstitutions']] = df[['aaSubstitutions']].fillna('')
    df[['aaDeletions']] = df[['aaDeletions']].fillna('')
    df[['aaInsertions']] = df[['aaInsertions']].fillna('')

    df['aa_mutations'] = df[['aaSubstitutions', 'aaDeletions', 'aaInsertions']].apply(lambda x: ','.join(x), axis=1)

    for index,row in df.iterrows():
        temp_mut = row['aa_mutations']
        temp_mut = temp_mut.split(',')
        aa_mutations.append(temp_mut)

    return aa_mutations


def training(model_name,max_number_strain):               # Training the Word2Vec Model
    cores=multiprocessing.cpu_count()
    w2v_model = Word2Vec(min_count = 1,
                        window=max_number_strain,     
                        vector_size=300,
                        sample=6e-5,
                        alpha=0.03,
                        min_alpha=0.0007,
                        workers=cores,
                        sg=1)
    w2v_model.build_vocab(aa_mutations)
    w2v_model.train(aa_mutations,total_examples=w2v_model.corpus_count,epochs=100)
    w2v_model.save(f"{model_name}_epoch100.model")
    w2v_model.wv.save_word2vec_format(f"word2vec_{model_name}_epoch100.txt")


def get_vector(epi_id):                                   
    mut_dict = dict(zip(df['epi_id'],df['aa_mutations']))
    docvec = []
    vectors = []
    muts = mut_dict[epi_id].split(',')
    for mut in muts:
        vectors.append(manu_aa_model.get_vector(mut))     # Using the trained Word2Vec model to calculate amino acid mutation embeddings
    vectors = np.asarray(vectors)
    avg_vec = vectors.mean(axis=0)                        # Calculating viral strain embeddings
    docvec.append(avg_vec)
    docvec = np.asarray(docvec)
    return docvec



