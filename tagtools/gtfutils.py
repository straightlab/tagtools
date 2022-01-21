import subprocess
import csv
import numpy as np
import pickle

def make_tx_dict_fromCSV(annot_file):
    annot={}
    with open(annot_file,'r') as f:
        myreader=csv.reader(f, delimiter='\t')
        for row in myreader:
            annot[row[0]]=row[1:6]
    return annot

def make_exons_dict_fromGFF(gff_file, nmax=0):
    m={}
    
    cmd = "cat "+gff_file+""" | awk 'BEGIN{OFS="\t"}((!/^#/) && $3=="transcript")"""
    cmd+="""{x=$9; split(x,xarr,\";\"); for (i=1; i<=length(xarr); i++)"""
    cmd+=""" {split(xarr[i],yy,"="); zz[yy[1]]=yy[2];}; print zz["ID"], $4, $5, $7, $1}'"""
    
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split("\t") 
        T_id=read_data[0] #transcript I
        
        if T_id in m:
            print("oops, already found")
#             m[T_id]+=[read_data[3],read_data[1],read_data[2]]
        else:
            nok+=1
            readstrand=int(1) if read_data[3]=="+" else int(-1)
            m[T_id]=[readstrand,int(read_data[1]),int(read_data[2]),read_data[4]]
            
    p.terminate()
    del p
    
    m2={}
    cmd = "cat "+gff_file+""" | awk 'BEGIN{OFS="\t"}((!/^#/) && $3=="exon")"""
    cmd+="""{x=$9; split(x,xarr,\";\"); for (i=1; i<=length(xarr); i++)"""
    cmd+=""" {split(xarr[i],yy,"="); zz[yy[1]]=yy[2];}; print zz["Parent"], $4, $5}'"""
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split("\t") 
        T_id=read_data[0] #transcript I
        
        vals=m2.pop(T_id,[])
        if len(vals)>0:
            
            m2[T_id]=np.vstack([vals,np.array([int(read_data[1]),int(read_data[2])])])
        else:
            nok+=1
            m2[T_id]=np.array([[int(read_data[1]),int(read_data[2])]])
    
    p.terminate()
    del p
    
    for k, v in m.items():
        vals=m2.pop(k,[])
        if len(vals)>0:
            if v[0]==1:
                a=np.argsort(vals[:,0])
                m2[k]=[v[3],np.insert(np.cumsum(vals[a,1]-vals[a,0]+1),0,0),vals[a,0]]
            else:
                a=np.argsort(-vals[:,1])
                m2[k]=[v[3],np.insert(np.cumsum(vals[a,1]-vals[a,0]+1),0,0),-vals[a,1]]
        else:
            m2[k]=[]
    
    return m2

def make_genebodies_dict_fromGFF(gff_file, nmax=0):
    m={}
    
    cmd = "cat "+gff_file+""" | awk 'BEGIN{OFS="\t"}((!/^#/) && $3=="gene")"""
    cmd+="""{x=$9; split(x,xarr,\";\"); for (i=1; i<=length(xarr); i++)"""
    cmd+=""" {split(xarr[i],yy,"="); zz[yy[1]]=yy[2];}; print zz["ID"], $4, $5, $7, $1}'"""
    
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split("\t") 
        T_id=read_data[0] #transcript I
        
        if T_id in m:
            print("oops, already found")
#             m[T_id]+=[read_data[3],read_data[1],read_data[2]]
        else:
            nok+=1
            readstrand=int(1) if read_data[3]=="+" else int(-1)
            m[T_id]=[readstrand,int(read_data[1]),int(read_data[2]),read_data[4]]
            
    p.terminate()
    del p
    
    m2={}
    cmd = "cat "+gff_file+""" | awk 'BEGIN{OFS="\t"}((!/^#/) && $3=="gene")"""
    cmd+="""{x=$9; split(x,xarr,\";\"); for (i=1; i<=length(xarr); i++)"""
    # cmd+=""" {split(xarr[i],yy,"="); zz[yy[1]]=yy[2];}; print zz["gene_id"], $4, $5}'"""
    # change on sept 16 2019 to avoid _PAR_Y genes to take over
    cmd+=""" {split(xarr[i],yy,"="); zz[yy[1]]=yy[2];}; print zz["ID"], $4, $5}'"""
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    nok=0
    for _, line in enumerate(p.stdout):
        if nmax>0 and nok>nmax:
            break
        read_data=line.strip().split("\t") 
        T_id=read_data[0] #transcript I
        
        vals=m2.pop(T_id,[])
        if len(vals)>0:
            
            m2[T_id]=np.vstack([vals,np.array([int(read_data[1]),int(read_data[2])])])
        else:
            nok+=1
            m2[T_id]=np.array([[int(read_data[1]),int(read_data[2])]])
    
    p.terminate()
    del p
    
    for k, v in m.items():
        vals=m2.pop(k,[])
        if len(vals)>0:
            if v[0]==1:
                a=np.argsort(vals[:,0])
                m2[k]=[v[3],np.insert(np.cumsum(vals[a,1]-vals[a,0]+1),0,0),vals[a,0]]
            else:
                a=np.argsort(-vals[:,1])
                m2[k]=[v[3],np.insert(np.cumsum(vals[a,1]-vals[a,0]+1),0,0),-vals[a,1]]
        else:
            m2[k]=[]
    
    return m2

def save_dict_toPickle(filename, d):
    with open(filename, 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_dict_fromPickle(filename):
    with open(filename, 'rb') as handle:
        d=pickle.load(handle)
        return d


def make_chr_length_dict_fromCSV(chr_file, delim='\t'):
    chrdict = {}
    with open(chr_file, 'r') as f:
        myreader = csv.reader(f, delimiter=delim)
        for _, row in enumerate(myreader):
            chrdict[row[0]] = int(row[1])
    return chrdict

def make_dict_fromCSV(annot_file, delim='\t', awk_filt=""):
    m={}
    cmd = "cat "+annot_file
    if len(awk_filt):
        cmd+=" | awk "+awk_filt
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True, start_new_session=True)
    for _, line in enumerate(p.stdout):
        myvals=line.strip().split(delim) 
        m[myvals[0]]=myvals[1:]

    p.terminate()
    del p
    return m
