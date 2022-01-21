import csv
import json
import pysam
import numpy as np
import pandas as pd
import random
from .bamutils import fetch_upto_next_ID
from . import gtfutils
DEFAULT_ANNOT=['?','?','?','?',int(0)]

def model(file_ambiv_gps, file_expression_model):
        #load ambivalence gps
    with open(file_ambiv_gps) as f:
        json_data=f.read()
    
    ambiv_gps_tx = json.loads(json_data)[0]
    ambiv_gps_tx_LUT = {v:k for k, v in ambiv_gps_tx.items()}
    
    #load expression model (hash table)
    if not file_expression_model is None:
        expression_model=_load_expression_model_sf(file_expression_model)

        #compute the probabilities
        ambivalence_proba, tx_LUT, tx_ix_dict = _compute_proba(ambiv_gps_tx_LUT, expression_model)
    else:
        ambivalence_proba, tx_LUT, tx_ix_dict = _compute_proba(ambiv_gps_tx_LUT, {})

    return ambivalence_proba, tx_LUT, tx_ix_dict
    # self.ambivalence_proba=ambivalence_proba
    # self.tx_LUT=tx_LUT
    # self.tx_ix_dict=tx_ix_dict

def _load_expression_model_sf(file, ix=3): #load a salmon file model
    with open(file,'r') as f:
        myreader=csv.reader(f, delimiter='\t')
        myreader.__next__() #get rid of header
        expression_model={row[0]:float(row[ix]) for row in myreader if float(row[ix])>0}
    return expression_model

def load_sf(file,ix): #load a salmon file model
    with open(file,'r') as f:
        myreader=csv.reader(f, delimiter='\t')
        myreader.__next__() #get rid of header
        expression_model={row[0]:float(row[ix]) for row in myreader if float(row[ix])>0}
    return expression_model



def _compute_proba(ambiv_gps_LUT, expression_model, tx_ix_dict=None): 
    # this assign probabibilities to all tx within ambivalence groups

    if tx_ix_dict is None:
        tx_ix_dict={}
        
    if len(tx_ix_dict)==0:
        i=-1
    else:
        i=max(v for _, v in tx_ix_dict.items())

    ambiv_proba=[[[],[]]]*(len(ambiv_gps_LUT)+1)
    for k, v in ambiv_gps_LUT.items():
        tx_all=v.split("|")
        tx_all_proba=[expression_model.get(tx,0) for tx in tx_all]
        tx_all_ix=[0]*len(tx_all_proba)
        for ix, tx in enumerate(tx_all):
            tx_ix=tx_ix_dict.get(tx,-1)
            if tx_ix<0:
                i+=1
                tx_ix_dict[tx]=i
                tx_all_ix[ix]=i
            else:
                tx_all_ix[ix]=tx_ix
        ambiv_proba[k]=[tx_all_proba,tx_all_ix]
    tx_LUT_dict={v:k for k,v in tx_ix_dict.items()} #{1:ENST*, 2:ENSTXXX,...}
    tx_LUT=[""]*(i+1)
    for k in range(i+1):
        tx_LUT[k]=tx_LUT_dict.get(k)

    return ambiv_proba, tx_LUT, tx_ix_dict


class Quantifier: #quantify .tag.bam file

    def __init__(self, file_ambiv_gps, file_expression_model=None): #agi: ambivalence group index
        ambivalence_proba, tx_LUT, tx_ix_dict=model(file_ambiv_gps, file_expression_model)
        self.ambivalence_proba=ambivalence_proba
        self.tx_LUT=tx_LUT
        self.tx_ix_dict=tx_ix_dict
        self.annot_dict=None
    
    def set_annotations_dictionary(self, annot_file):
        # Load annotation dictionnary
        if isinstance(annot_file,str):
            print("Loading annotations dictionary from %s"%annot_file)
            self.annot_dict=gtfutils.make_tx_dict_fromCSV(annot_file)
        else:
            print('Annotation dictionnary passed directly as argument')
            self.annot_dict=annot_file
    

    def quant(self, qbam_in, todf=True, simplify=True, drop_zeros=True, allow_no_prior=True):
        # dumps a quantification table
        #prepares the output table for know tx, and an empty dict for unknown stuff
        noprior=1.0 if allow_no_prior else 0
        q=np.zeros((len(self.tx_ix_dict),10))
        q_unknown={}
        q_unknown2={}
        n_ambiv_with_prior=0
        n_ambiv_without_prior=0
        n_unique_with_prior=0
        n_unique_new=0

        reader = pysam.AlignmentFile(qbam_in)
        bamiter = reader.fetch(until_eof=True)

        try:
            r = bamiter.__next__()
        except StopIteration: 
            r = None
        while not(r==None):
            ambiv_ix = r.get_tag("ai") if r.has_tag("ai") else 0
            ambiv_ix_abs=abs(ambiv_ix)
            ambiv_ix_gene = r.get_tag("aI") if r.has_tag("aI") else 0
            ambiv_gene_offset = 0 if abs(ambiv_ix_gene)>1 else 3 #3 if NOT ambiv at gene level
            offset=0 if ambiv_ix>0 else 5 #if negative ix, means there is other alignment
            # q cols: 0 is ambiv_tx with prior, 1 is ambiv_tx no prior, 2 is no ambiv, 
            # 3 is ambiv at tx but not at gene, 4 same no prior

            if ambiv_ix_abs>1: #this is an ambivalence group
                Es, tx_ixs =self.ambivalence_proba[abs(ambiv_ix)] #expression values 
                Etotal=sum(Es)
                if Etotal>0:
                    n_ambiv_with_prior+=1
                    for i, tx_ix in enumerate(tx_ixs):
                        q[tx_ix,0+ambiv_gene_offset+offset]+=Es[i]/Etotal
                else:
                    n_ambiv_without_prior+=1
                    for i, tx_ix in enumerate(tx_ixs):
                        q[tx_ix,1+ambiv_gene_offset+offset]+=1.0/len(tx_ixs)
            elif r.has_tag("ar"):
                
                tx_name=r.get_tag("ar")
                tx_ix=self.tx_ix_dict.get(tx_name,-1)
                if tx_ix>-1: # this is a known one
                    n_unique_with_prior+=1
                    q[tx_ix,2+offset]+=1
                else:
                    n_unique_new+=1
                    if offset==0:
                        if tx_name in q_unknown: 
                            q_unknown[tx_name]+=1
                        else:
                            q_unknown[tx_name]=1
                    else:
                        if tx_name in q_unknown2: 
                            q_unknown2[tx_name]+=1
                        else:
                            q_unknown2[tx_name]=1
            
            try:
                r = bamiter.__next__()
            except StopIteration: 
                r = None

        reader.close()
        counts=(n_ambiv_with_prior,n_ambiv_without_prior,n_unique_with_prior,n_unique_new)
        print("Found: %g ambiv_with_prior, %g ambiv_without_prior, %g unique known, %g unknown"%counts)

        #Now consolidate

        #create LUT dict for the new tx
        tx_new_LUT_dict={i:k for i,(k,_) in enumerate(q_unknown.items())} #{1:ENST*, 2:ENSTXXX,...}
        imax=len(tx_new_LUT_dict)
        #add in the new tx that belong to the ai=-1 gp
        for k, _ in q_unknown2.items():
            if not k in q_unknown:
                tx_new_LUT_dict[imax]=k
                imax+=1
        
        #transform the LUT dict into a list
        tx_new_LUT=[""]*(len(tx_new_LUT_dict))
        for k in range(len(tx_new_LUT_dict)):
            tx_new_LUT[k]=tx_new_LUT_dict.get(k)
        
        #transform the q_unknown dictionnary into a list
        q_new=np.zeros((len(tx_new_LUT_dict),2))
        for k, tx in enumerate(tx_new_LUT):
            counts=q_unknown.get(tx,0)
            counts2=q_unknown2.get(tx,0)
            q_new[k,0]=counts
            q_new[k,1]=counts2

        
        if not todf:
            return q, q_new, tx_new_LUT, counts
        else:
    
            #ambiv_tx is ambiv at tx level (not at gene)
            #ambiv_gene is ambiv also at gene level
            #no prior means no_prior proba for this ambivalence group
            #2 means there is alignment to other areas

            # df_known=pd.DataFrame(q,columns=['ambiv','ambiv_noprior','unique','ambiv2','ambiv_noprior2','unique2'], index=self.tx_LUT)
            # df_unknown=pd.DataFrame(q_new,columns=['unique','unique2'], index=tx_new_LUT)
            # df_all=pd.concat([df_known,df_unknown], axis=0)
            # df_all.index.names=['Name']
            # if simplify:
            #     df=pd.concat(
            #         [df_all.sum(axis=1),df_all["ambiv"]+df_all["ambiv2"],df_all["ambiv2"]+df_all["unique2"]],axis=1).drop(["ENST*"])
            #     df.columns=["N","N_ambiv","N_hasunannot"]
            #     if drop_zeros:
            #         df.fillna(0, inplace=True)
            #         df.drop(df[(df.N ==0)].index, inplace=True)
            #     return df
            # else:
            #     return df_all
            df_known=pd.DataFrame(q,columns=['ambiv_gene','ambiv_gene_noprior','unique','ambiv_tx','ambiv_tx_noprior','ambiv_gene2','ambiv_gene_noprior2','unique2','ambiv_tx2','ambiv_tx_noprior2'], index=self.tx_LUT)
            df_unknown=pd.DataFrame(q_new,columns=['unique','unique2'], index=tx_new_LUT)
            df_all=pd.concat([df_known,df_unknown], axis=0)
            df_all.index.names=['Name']
            if simplify:
                df=pd.concat(
                    [df_all.sum(axis=1),
                    df_all["ambiv_tx"]+df_all["ambiv_tx2"]+noprior*(df_all["ambiv_tx_noprior"]+df_all["ambiv_tx_noprior2"]),
                    df_all["ambiv_gene"]+df_all["ambiv_gene2"]+noprior*(df_all["ambiv_gene_noprior"]+df_all["ambiv_gene_noprior2"]),
                    df_all["unique2"]+df_all["ambiv_tx2"]+df_all["ambiv_gene2"]+noprior*(df_all["ambiv_tx_noprior2"]+df_all["ambiv_gene_noprior2"]),
                    df_all["ambiv_tx2"]+noprior*(df_all["ambiv_tx_noprior2"]),
                    df_all["ambiv_gene2"]+noprior*(df_all["ambiv_gene_noprior2"])],axis=1).drop(["ENST*"])

                df.columns=["N","N_ambiv_tx", "N_ambiv_gene", "N_hasunannot","N_ambiv_tx_hasunannot", "N_ambiv_gene_hasunannot"] #"N_ambiv_tx" ambig at tx level (not at gene),  "N_ambiv_gene" (ambiv at gene) #in the simplification, allow or not no prior
                if drop_zeros:
                    df.fillna(0, inplace=True)
                    df.drop(df[(df.N ==0)].index, inplace=True)
                return df
            else:
                return df_all
    
    # def resolve():
    #     #resolves the ambivalence classes using either probabilisitic or top mode, this spits out another bam file

    def _resolve_bam_stream(self,bamiter, out):
        #go through file
        try:
            r = bamiter.__next__()
        except StopIteration: 
            r = None

        while not(r==None):
            ambiv_ix = r.get_tag("ai") if r.has_tag("ai") else 0
            ambiv_ix_abs=abs(ambiv_ix)
            # ambig=1 if ambiv_ix>0 else -1
            records, r = fetch_upto_next_ID(bamiter,r)

            if ambiv_ix_abs>1: #this is an ambivalence group
                Es, tx_ixs =self.ambivalence_proba[abs(ambiv_ix)] #expression values 
                Etotal=sum(Es)
                i=list(range(len(Es)))
                if Etotal>0: #it has prior proba
                    i_sampled=random.choices(i,Es)[0]
                    rec_sampled_ix=-1
                    for ii, rec in enumerate(records):
                        ambiv_ix_within= list(rec.get_tag("ay")) if rec.has_tag("ay") else None
                        if ambiv_ix_within is None:
                            ambiv_ix_within=[rec.get_tag("ax")] if rec.has_tag("ax") else []
                        # else:
                            # print(ambiv_ix_within)
                        if i_sampled in ambiv_ix_within:
                            rec_sampled_ix=ii
                            break
                    if rec_sampled_ix>-1:
                        
                        records[rec_sampled_ix].set_tag("aq",Es[i_sampled]/Etotal, value_type='f')
                        #update the tags to reflect possible change in representative annotation
                        if not self.annot_dict is None:
                            newtx=self.tx_LUT[tx_ixs[i_sampled]]
                            annot_desc=self.annot_dict.get(newtx, DEFAULT_ANNOT)
                            records[rec_sampled_ix].set_tag('an',annot_desc[1],value_type='Z')
                            records[rec_sampled_ix].set_tag('at',annot_desc[2],value_type='Z')
                            #records[rec_sampled_ix].set_tag('as',annot_desc[3],value_type='A')
                            records[rec_sampled_ix].set_tag('ag',annot_desc[0],value_type='Z')
                        
                        out.write(records[rec_sampled_ix])
                    else:
                        print("Warning: not able to pick tx for sampling")

                
                else: #no prior proba
                    i_sampled=random.choice(i)
                    rec_sampled_ix=-1
                    for ii, rec in enumerate(records):
                        ambiv_ix_within= list(rec.get_tag("ay")) if rec.has_tag("ay") else None
                        if ambiv_ix_within is None:
                            ambiv_ix_within=[rec.get_tag("ax")] if rec.has_tag("ax") else []
                        # else:
                            # print(ambiv_ix_within)
                        if i_sampled in ambiv_ix_within:
                            rec_sampled_ix=ii
                            break
                    if rec_sampled_ix>-1:
                       
                        records[rec_sampled_ix].set_tag("aq",-1.0/len(i), value_type='f')
                        out.write(records[rec_sampled_ix])
                    else:
                        print("Warning: not able to pick tx for sampling")                    
                   
            else: 
                if len(records)>1:
                    print("Warning: record with multiple alignments but no ambivalence group")
                    records[0].set_tag("aq",0.0, value_type='f')
                out.write(records[0])
                
    def resolve(self, inbam, outbam):
        bamreader = pysam.AlignmentFile(inbam)
        bamiter = bamreader.fetch(until_eof=True)
   
        if not outbam is None:
            outmode="wb" if outbam.endswith(".bam") else "w"
            out=pysam.AlignmentFile(outbam, outmode, template=bamreader)
        else:
            out=pysam.AlignmentFile("/dev/stdout", "w", template=bamreader)

        self._resolve_bam_stream(bamiter,out)

        bamreader.close()
        out.close()



    

    