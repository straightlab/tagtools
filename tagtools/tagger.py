#!/usr/bin/env python
import pysam
# from . import bamutils.*
import numpy as np
import sys
from . import gtfutils
import pickle
import json
import bisect
import array
DEFAULT_ANNOT=['*','*','*','*',int(0)]


def fetch_upto_next_ID(bam_reader,r): # r is the current record
    records_list=[r] #the records with the same readID
#     print(records_list[0].to_string())
    id_wanted=r.query_name
    try:
        record_current=bam_reader.__next__()
        id_current=record_current.query_name
        while id_current==id_wanted:

            records_list+=[record_current]
            record_current=bam_reader.__next__()
            id_current=record_current.query_name
        return records_list, record_current
    except StopIteration:
        return records_list, None

def fetch_ID(bam_reader,r,id_wanted):
    while ((not r is None) and r.query_name<id_wanted):
        try:
            r=bam_reader.__next__()
        except StopIteration:
            return [], None
        
    if not r is None and r.query_name==id_wanted:
        return fetch_upto_next_ID(bam_reader,r)
    else:
        return [], r


def _lt_natural(x,y, naturalSort):
    if naturalSort:
        return ([ int(c) if c.isdigit() else c.lower() for c in x.split(':') ] < [ int(c) if c.isdigit() else c.lower() for c in y.split(':') ] )
    else:
        return x<y

def _convert_pos_tx_to_genomic(exon_map,pos,s): 
    #pass reference end for negative align, reference_start +1 for positive
    ix=np.searchsorted(exon_map[1],pos)
    if ix>0 and ix<len(exon_map[1]):
        return s*(exon_map[2][ix-1]+pos-exon_map[1][ix-1]-1)
    else:
        return 0


def _genomic_position_from_tx_alignment(r, x):
    outpos=-1
    outchr="*"
    outpos=0
    outpos2=0
    if len(x)>0:
        outchr=x[0]
        
        if x[2][0]<0:
            pos=r.reference_end
            ix=np.searchsorted(x[1],pos)
            if ix>0 and ix<len(x[1]):
                outpos=-(x[2][ix-1]+pos-x[1][ix-1]-1)

            pos2=r.reference_start+1
            ix2=np.searchsorted(x[1],pos2)
            if ix2>0 and ix2<len(x[1]):
                outpos2=-(x[2][ix2-1]+pos2-x[1][ix2-1]-1)

            l=outpos2-outpos+1
            
        else:
            pos=r.reference_start+1
            ix=np.searchsorted(x[1],pos)
            if ix>0 and ix<len(x[1]):
                outpos=(x[2][ix-1]+pos-x[1][ix-1]-1)

            pos2=r.reference_end
            ix2=np.searchsorted(x[1],pos2)
            if ix2>0 and ix2<len(x[1]):
                outpos2=(x[2][ix2-1]+pos2-x[1][ix2-1]-1)

            l=outpos2-outpos+1
    
    return "_".join([outchr,str(outpos),str(l)])

def _pairup_genomic_and_tx_alignments(records_ref, records_annot, exons_dict, isunstranded, isrevstrand):
    # there is a single primary alignment
    nannot=len(records_annot)
    pairs_props=np.zeros((nannot,4),int)

    #CASE where this is an unmapped read, which should be present in both the genome and tx bam, single record
    if nannot==1 and (records_annot[0].is_unmapped): 
        if len(records_ref)>1: 
            print('ReadID %s has multiple unaligned entries, something is wrong '%records_annot[0].query_name)
        return [], 0, 0
    
    #CASE this is mapped read
    else:
    
        # ref_dict={(r.reference_name+"_"+str(r.reference_start+1-r.query_alignment_start)):ix for ix, r in enumerate(records_ref)}
        # changed on sep 16 2019
        ref_dict={(r.reference_name+"_"+str(r.reference_start+1)+"_"+str(r.reference_end-r.reference_start)):ix for ix, r in enumerate(records_ref)}

        ixs=[ref_dict.get(_genomic_position_from_tx_alignment(r_annot, exons_dict.get(r_annot.reference_name,[])),-1) 
             for r_annot in records_annot]
        
        pairs_props[:,0]=ixs
        primary=-1 # the primary 
        # nonsupp=-1
        sortixs=np.argsort(pairs_props[:,0])
        hasannot=[False]*len(records_ref)
        hasannot_strandok=[False]*len(records_ref)

        for _, i in enumerate(sortixs):
            
            if pairs_props[i,0]<0: #this should not occur
                print('Found annotation without genomic counterpart: %s _ %s'%(
                    records_annot[i].query_name, _genomic_position_from_tx_alignment(
                    records_annot[i], exons_dict.get(records_annot[i].reference_name,[]))))
                print(ref_dict)
            else:
                hasannot[pairs_props[i,0]]=True
                if (isunstranded or records_annot[i].is_reverse==isrevstrand):
                    hasannot_strandok[pairs_props[i,0]]=True
                    pairs_props[i,3]=1
                    if (primary<0 or primary==pairs_props[i,0]): # this makes the first genomic alignment wich has an annotation on a valid strand the primary 
                        # alignment (non secondary). Here, the term primary refers to the genomic value of the alignment. All the annotations with the same genomic alignment will also by primary, but only one will be represenstative, the others will be supplementatry
                        pairs_props[i,1]=1
                        primary=pairs_props[i,0]

        reprannot=-1 #representative annotation
        for _, i in enumerate(sortixs):
            refix=pairs_props[i,0]
            if reprannot<refix and ((hasannot[refix]==False) or (hasannot_strandok[refix]==False) or (pairs_props[i,3]==1)):
                pairs_props[i,2]=1 #make it representative (non supplementary). Only ONE annotation per genomic alignment is representative
                reprannot=refix


                    
                    # hasannot_strandwrong[pairs_props[i,0]]=True
        n_noannot=len(records_ref)-np.sum(hasannot_strandok) #n_noannot includes genomic alignment with no annotation at all or with annotations on the wrong strand
        n_hasannot_wrongstrand=np.sum(hasannot)-np.sum(hasannot_strandok)

        #
        # n_hasannot_okstrand=np.sum(hasannot_strandok)
        # if n_hasannot_wrongstrand>0:
        #     print('Here')
        #     print(pairs_props)
        #     print('There')
        return pairs_props, n_noannot, n_hasannot_wrongstrand



def _compute_ambivalence_groups(records_annot, pairs_props, ambivalence_dict, ixmax, annot_dict, ambivalence_dict_GENE, ixmax_GENE):
    
    ambivalence_group_ix=1 #one means no ambivalence
    ambivalence_group_ix_GENE=1 #one means no ambivalence

    n_strandok=np.sum(pairs_props[:,3]==1) #number of GENOMIC alignments that have at least one annotation on the proper strand
    
    annots_ag_ix=-1*np.ones((pairs_props.shape[0],2),int) #index within the ambivalence group

    if n_strandok>1 and (not ambivalence_dict is None): #we are dealing with an ambivalent annotation
        # if np.all(pairs_props[:,3]==0):
        #     print('Computing an equivalence group for geneid %s'%records_annot[0].query_name)
        #     print(pairs_props)
        #     print("l=%g"%len(records_annot))
        ambivalence_list=list(set([records_annot[i_annot].reference_name for i_annot, ix in enumerate(pairs_props) if (ix[3]==1)]))
        ambivalence_list.sort()
        nag=len(ambivalence_list)
       
        for i, r in enumerate(records_annot):
            if pairs_props[i,3]==1:
                pos = bisect.bisect_left(ambivalence_list, r.reference_name)
                if pos != nag and ambivalence_list[pos] == r.reference_name:
                    annots_ag_ix[i,0]=pos

        ambivalence_group_name="|".join(ambivalence_list)
        ambivalence_group_ix=ambivalence_dict.get(ambivalence_group_name,0)
        if ambivalence_group_ix==0:
            ixmax+=1
            ambivalence_dict[ambivalence_group_name]=ixmax
            ambivalence_group_ix=ixmax

        if not ambivalence_dict_GENE is None:
            genes_list=[annot_dict.get(tx,'*')[0] for tx in ambivalence_list]
            ambivalence_list_GENE=list(set(genes_list))
            if len(ambivalence_list_GENE)>1:
                ambivalence_list_GENE.sort()
                
                nag_GENE=len(ambivalence_list_GENE)
       
                for i in range(len(records_annot)):
                    if annots_ag_ix[i,0]>-1:
                        pos = bisect.bisect_left(ambivalence_list_GENE, genes_list[annots_ag_ix[i,0]])
                        if pos != nag_GENE and ambivalence_list_GENE[pos] == genes_list[annots_ag_ix[i,0]]:
                            annots_ag_ix[i,1]=pos

                ambivalence_group_name_GENE="|".join(ambivalence_list_GENE)
                ambivalence_group_ix_GENE=ambivalence_dict_GENE.get(ambivalence_group_name_GENE,0)
                if ambivalence_group_ix_GENE==0:
                    ixmax_GENE+=1
                    ambivalence_dict_GENE[ambivalence_group_name_GENE]=ixmax_GENE
                    ambivalence_group_ix_GENE=ixmax_GENE
        
    return ambivalence_group_ix, ixmax, ambivalence_group_ix_GENE, ixmax_GENE, annots_ag_ix

def _crossannotate_alignments_and_write(records_ref, records_annot, exons_dict, annot_dict, ambivalence_dict, ixmax, ambivalence_dict_GENE, ixmax_GENE, bam_out, bam_out_noannot, isunstranded, isrevstrand, updateflags, keep_noannot, collapse_annotations):
    
    pairs_props, n_noannot, n_hasannot_wrongstrand=_pairup_genomic_and_tx_alignments(records_ref, records_annot, exons_dict, isunstranded, isrevstrand)
    # n_noannot is the number of genomic matches with no alignment on the proper strand
    # CASE this is an unmapped read, don't change anything
    if len(pairs_props)==0: 
        n_with_annot_strandok=0
        n_with_annot_strandwrong=0
        n_nomatch=0
        for rec_ref in records_ref: # there should be only one read in records_ref
            if keep_noannot:
                bam_out_noannot.write(rec_ref)
            n_nomatch+=1
        if n_nomatch>1:
            print('Found read with multiple no match, something is wrong')
       
    
    # CASE this is a mapped read
    else:    

        n_with_annot_strandok=np.sum(pairs_props[:,3]) #this is the total number of annot on the proper strand
        n_with_annot_strandwrong=pairs_props.shape[0]-n_with_annot_strandok #this is the number of annotations on the wrong strand
        n_nomatch=0 
        
        ambivalence_group_ix, ixmax, ambivalence_group_ix_GENE, ixmax_GENE, annots_ag_ix = _compute_ambivalence_groups(records_annot, pairs_props, ambivalence_dict, ixmax, annot_dict, ambivalence_dict_GENE, ixmax_GENE)
        
        # we multiply the ambivalence groups by -1 if there are genomic alignmnets with no annotation on the correct strand
        if n_noannot>0:
            ambivalence_group_ix_GENE=-ambivalence_group_ix_GENE
            ambivalence_group_ix=-ambivalence_group_ix
        
        # annotate the reads
        for i, rec_ref in enumerate(records_ref):
            
            flag_val=rec_ref.flag
            
            if collapse_annotations==True:
                # added & (pairs_props[:,3]==1) on nov 26 to make sure we only report annotations on the right strand
                ixs_annot=np.argwhere((pairs_props[:,0]==i) & (pairs_props[:,2]==1) & (pairs_props[:,3]==1) ).flatten() #only keep representative annotation (even a genomic alignment without a sense annot has a representative annot "none" or the wrong strand)
                # print(len(np.argwhere(pairs_props[:,0]==i).flatten()),len(ixs_annot))
                ixs_annot_all=np.argwhere(pairs_props[:,0]==i).flatten()
                annots_ag_ix_list=set(annots_ag_ix[ixs_annot_all,0]) #for this genomic alignment, list all the
                annots_ag_ix_list.discard(-1)
                annots_ag_ix_GENE_list=set(annots_ag_ix[ixs_annot_all,1])
                annots_ag_ix_GENE_list.discard(-1)
            else:
                ixs_annot=np.argwhere(pairs_props[:,0]==i).flatten()
                annots_ag_ix_list=[] #no need to keep track for uncompressed mode

            # CASE this genomic alignmnent has no annotation at all on the right strand (or this is simply an unmapped)
            if (len(ixs_annot)==0) and keep_noannot and (n_with_annot_strandok==0) and (i==0): # only report this one if there are no alignment on the right strand, and only keep one of the alignment
                rec_ref.set_tag('gS',len(records_ref)-n_noannot,value_type='i') # genomic SENSE of genomic alignments with annotations on the proper strand (and possibly on wrong strand)
                rec_ref.set_tag('gN',n_noannot,value_type='i') # number of genomic alignments with annotations on the wrong strand
                rec_ref.set_tag('aS',n_with_annot_strandok,value_type='i') # total number of possible annot with sense orientation (redundant with ag but that's ok)
                rec_ref.set_tag('ar','*',value_type='Z')
                rec_ref.set_tag('ap',int(0),value_type='i')
                rec_ref.set_tag('ai',ambivalence_group_ix,value_type='i') #"ambivalence index"
                rec_ref.set_tag('aI',ambivalence_group_ix_GENE,value_type='i') #"gene ambivalence index"

                if updateflags:
                    rec_ref.flag = flag_val | (1<<8) #set bit 256=secondary alignment
                    rec_ref.flag = rec_ref.flag | (1<<11) #set it as a NON representative alignment (representative alignment)
                
                
                bam_out_noannot.write(rec_ref)
           
           # CASE this genomic alignment has at least one annotation on the proper strand (and thus one of them is representative, preferentially a sense annotation if there is one, but an antisense one ow)
            else:
                for _, annot_id in enumerate(ixs_annot):
                    
                    rec_ref.set_tag('gS',len(records_ref)-n_noannot,value_type='i') # genomic SENSE of genomic alignments with annotations on the proper strand (and possibly on wrong strand)
                    rec_ref.set_tag('gN',n_noannot,value_type='i') # number of genomic alignments with annotations on the wrong strand or no annotations
                    rec_ref.set_tag('aS',n_with_annot_strandok,value_type='i') # total number of possible annot with sense orientation (redundant with ag but that's ok)
                    rec_ref.set_tag('ar',records_annot[annot_id].reference_name,value_type='Z')
                    rec_ref.set_tag('ap',int(records_annot[annot_id].reference_start),value_type='i')

                    rec_ref.set_tag('ai',ambivalence_group_ix,value_type='i') #"ambivalence index"
                    rec_ref.set_tag('aI',ambivalence_group_ix_GENE,value_type='i') #"gene ambivalence index"
                    if abs(ambivalence_group_ix)>1:
                        rec_ref.set_tag('ax',annots_ag_ix[annot_id,0],value_type='i') #within the ambivalence, what is the index
                    if abs(ambivalence_group_ix_GENE)>1:
                        rec_ref.set_tag('aX',annots_ag_ix[annot_id,1],value_type='i') #within the gene ambivalence, what is the index
                    if len(annots_ag_ix_list)>1:
                        # print(annots_ag_ix_list)
                        rec_ref.set_tag('ay',array.array('H',annots_ag_ix_list))
                    if len(annots_ag_ix_GENE_list)>1:
                        rec_ref.set_tag('aY',array.array('H',annots_ag_ix_GENE_list))
                    if not(annot_dict==None):
                        annot_desc=annot_dict.get(records_annot[annot_id].reference_name,DEFAULT_ANNOT)
                        rec_ref.set_tag('an',annot_desc[1],value_type='Z')
                        rec_ref.set_tag('at',annot_desc[2],value_type='Z')
                        #rec_ref.set_tag('as',annot_desc[3],value_type='A')
                        rec_ref.set_tag('ag',annot_desc[0],value_type='Z') #ENSG
                        # rec_ref.set_tag('al',int(annot_desc[4]),value_type='i')
                    
                    if updateflags:
                        rec_ref.flag=flag_val # reset to initial value
                        rec_ref.flag = (flag_val & ~(1<<8)) if pairs_props[annot_id,1]==1 else (flag_val | (1<<8)) # unset 256 if primary, set if secondary 
                        rec_ref.flag = (rec_ref.flag &~ (1<<11)) if pairs_props[annot_id,2]==1 else (rec_ref.flag | (1<<11)) # unset 2048 if representative else set (supplementary)

                    if pairs_props[annot_id,3]==1:
                        bam_out.write(rec_ref)
                    #changed nov 26 2019, now we don't keep the no annotation if there is
                    # elif keep_noannot and (n_with_annot_strandok==0):
                    #     bam_out_noannot.write(rec_ref)

    return n_with_annot_strandok, n_with_annot_strandwrong, n_noannot, n_hasannot_wrongstrand, n_nomatch, ixmax, ixmax_GENE


def _annotate_BAM_stream(bamiter_ref,bamiter_annot,exons_dict, annot_dict, ambivalence_dict, ixmax, ambivalence_dict_GENE, ixmax_GENE, bam_out, bam_out_noannot, isunstranded, 
                              isrevstrand, updateflags, keep_noannot=True, collapse_annotations=True, maxd=200, maxd_genomic=20, nmax=0):
    
    fulloutput=True
    if nmax>0:
        fulloutput=False
    
    #annotation level stats
    stats_annot_SENSE=[int(0) for i in range(maxd+1)] #total number of annotation on ok strand
    stats_annot_ANTISENSE=[int(0) for i in range(maxd+1)] #total number of annotation on wrong strand
    stats_annot_ANY=[int(0) for i in range(maxd+1)]
    annot_SENSE=0
    annot_ANTISENSE=0

    #genomic mapping level stats 
    stats_gmappings=[int(0) for i in range(maxd_genomic+1)] #number of genomic alignments
    stats_gmappings_annotSENSE=[int(0) for i in range(maxd_genomic+1)] #number of genomic alignments with annoation on the correct strand
    stats_gmappings_annotANTISENSE=[int(0) for i in range(maxd_genomic+1)] #number of genomic alignments with annotaion on the wrong strand only 
    stats_gmappings_annotNONE=[int(0) for i in range(maxd_genomic+1)]#number of genomic alignments without annotation at all
    stats_gmappings_unanottatedONLY=[int(0) for i in range(maxd_genomic+1)]#number of genomic alignmnets with NO annotation of the correct strat
    gmappings=0
    gmappings_annotSENSE=0
    gmappings_annotANTISENSE=0
    gmappings_annotNONE=0
    gmappings_unanottatedONLY=0

    
    #read level stats
    tot_reads=0
    tot_noref=0 #should stay at 0
    reads_mapped_SENSEannot_unique=0 #unique alignmnet wiht sense annot
    reads_mapped_SENSEannot_multiple=0 #mulitple alignments with sense annot
    reads_mapped_ANTISENSEannot=0 #all alignments have ANTISENSE MAPPING OR NONE
    reads_mapped_NOannot=0 # all alignments have No annotation
    reads_unmapped=0
    # load the first read to initialize
    try:
        r_ref = bamiter_ref.__next__()
    except StopIteration: 
        r_ref = None
    try:
        r_annot = bamiter_annot.__next__()
    except StopIteration: 
        r_annot= None

    while not(r_ref==None) and (fulloutput or tot_reads<nmax):
        if r_annot==None or (r_ref.query_name<r_annot.query_name): # no more annotation data, or ref stream is behind, so we write the reference alone. The ambiguity is 0 (not labeled)

            records_ref, r_ref = fetch_upto_next_ID(bamiter_ref,r_ref)
            n_noannot=len(records_ref) #these have no annotation at all
            n_unmapped=0

            for i, rec_ref in enumerate(records_ref):
                if rec_ref.is_unmapped:
                    n_unmapped+=1
                    n_noannot=0
                    if keep_noannot and (bam_out_noannot is None) and (i==0):
                        bam_out.write(rec_ref)
                    elif keep_noannot and (i==0):
                        bam_out_noannot.write(rec_ref)
                    
                else:
                    
                    if updateflags:
                        rec_ref.flag |= (1<<8) #set bit 256=secondary alignment

                    rec_ref.set_tag('gS',0,value_type='i') # number of genomic alignments with annotations on the proper strand (and possibly on wrong strand)
                    rec_ref.set_tag('gN',0,value_type='i') # number of genomic alignments with no annotation at all or with annotations on the wrong strand
                    #the number of alignments that hve NO annotations at all is thus NH-nP-nM
            
                    # NUMBER OF CORRECT ANNOTATIONS
                    rec_ref.set_tag('aS',0,value_type='i') # total number of possible annot with correct orientation (redundant with ag but that's ok)
                    rec_ref.set_tag('ar','*',value_type='Z')
                    rec_ref.set_tag('ap',int(0),value_type='i')
                    rec_ref.set_tag('ai',0,value_type='i') #"ambivalence index"
                    rec_ref.set_tag('aI',0,value_type='i') #"gene ambivalence index"

                    if keep_noannot and (bam_out_noannot is None) and (i==0):
                        bam_out.write(rec_ref)
                    elif keep_noannot and (i==0):
                        bam_out_noannot.write(rec_ref)
            
            tot_reads+=1

            gmappings_unanottatedONLY+=n_noannot #this doesn't get incremented if this is a no match
            gmappings+=n_noannot
            stats_gmappings[min(n_noannot,maxd_genomic)]+=1
            if n_noannot>0: #only update when this is an actual matched read
                stats_gmappings_unanottatedONLY[min(n_noannot,maxd_genomic)]+=1
                reads_mapped_NOannot+=1 
            else:
                reads_unmapped+=1

    
        elif r_ref.query_name==r_annot.query_name: #this is a match SO at this point, the ref and annot streams are pointing to unanalzyed records with matching IDs,

            records_annot, r_annot = fetch_upto_next_ID(bamiter_annot,r_annot)
            records_ref, r_ref = fetch_upto_next_ID(bamiter_ref,r_ref)
            
            n_annot_strandok, n_annot_strandwrong, n_noannot, n_hasannot_wrongstrandONLY, n_unmapped, ixmax, ixmax_GENE = _crossannotate_alignments_and_write(records_ref, records_annot, exons_dict, annot_dict, ambivalence_dict, ixmax, ambivalence_dict_GENE, ixmax_GENE, bam_out, bam_out_noannot, isunstranded, isrevstrand, updateflags, keep_noannot, collapse_annotations)
            
            # update read level stats
            n_ref=len(records_ref)
            tot_reads+=1

            #update genomic alignment level stats
            stats_gmappings[min(n_ref-n_unmapped,maxd_genomic)]+=1 #number of matches (number of genomic alignmnents, no including nomatch)
            gmappings+=(n_ref-n_unmapped)
          
            n_hasannot_strandok=n_ref-n_unmapped-n_noannot
            stats_gmappings_annotSENSE[min(n_hasannot_strandok,maxd_genomic)]+=1 #number of genomic alignments with at least one annotation on the ok strand
            gmappings_annotSENSE+=n_hasannot_strandok

            stats_gmappings_annotANTISENSE[min(n_hasannot_wrongstrandONLY,maxd_genomic)]+=1 #number of genomic alignments with annotations only on opposite strand
            gmappings_annotANTISENSE+=n_hasannot_wrongstrandONLY

            n_gmappings_annotNONE=n_noannot-n_hasannot_wrongstrandONLY
            stats_gmappings_annotNONE[min(n_gmappings_annotNONE,maxd_genomic)]+=1 #number of genomic alignments with no annotation (but some with read have)
            gmappings_annotNONE+=(n_gmappings_annotNONE)

            #annotation level stats
            stats_annot_ANY[min(n_annot_strandok+n_annot_strandwrong,maxd)]+=1 #number of annotation (plus or minus)
            
            stats_annot_SENSE[min(n_annot_strandok,maxd)]+=1
            annot_SENSE+=n_annot_strandok
        
            stats_annot_ANTISENSE[min(n_annot_strandwrong,maxd)]+=1
            annot_ANTISENSE+=n_annot_strandwrong

            #read level stats
            if n_hasannot_strandok==1:
                reads_mapped_SENSEannot_unique+=1
            elif n_hasannot_strandok>1:
                reads_mapped_SENSEannot_multiple+=1
            elif n_hasannot_wrongstrandONLY>0:
                reads_mapped_ANTISENSEannot+=1
            elif n_unmapped==0:
                reads_mapped_NOannot+=1
            else:
                reads_unmapped+=1
        else: #no match an the annotation stream is behind, so just move the annotation cursor up, this should not happen
            tot_noref+=1
            try:
                r_annot = bamiter_annot.__next__()
            except StopIteration:
                r_annot=None

    d_tot={'reads':tot_reads, 'reads_mapped_SENSEannot_unique': reads_mapped_SENSEannot_unique,'reads_mapped_SENSEannot_multiple': reads_mapped_SENSEannot_multiple,
    'reads_mapped_ANTISENSEannot':reads_mapped_ANTISENSEannot, 'reads_mapped_NOannot':reads_mapped_NOannot, 'reads_unmapped':reads_unmapped, 'gmappings':gmappings, 'noref':tot_noref, 
       'gmappings_annotSENSE':gmappings_annotSENSE, 'gmappings_annotANTISENSE':gmappings_annotANTISENSE, 'gmappings_annotNONE':gmappings_annotNONE,
       'gmappings_unanottatedONLY':gmappings_unanottatedONLY, 'annot_SENSE':annot_SENSE, 'annot_ANTISENSE':annot_ANTISENSE}
    
    
    d_stats={'stats_gmappings':stats_gmappings,'stats_gmappings_annotSENSE':stats_gmappings_annotSENSE, 'stats_gmappings_annotANTISENSE':stats_gmappings_annotANTISENSE,
       'stats_gmappings_annotNONE':stats_gmappings_annotNONE, 'stats_gmappings_unanottatedONLY':stats_gmappings_unanottatedONLY, 'stats_annot_ANY':stats_annot_ANY,  'stats_annot_SENSE': stats_annot_SENSE,
       'stats_annot_ANTISENSE':stats_annot_ANTISENSE}

    return d_tot, d_stats

#exons_dict, annot_dict, bam_out, isunstranded, isrevstrand
def annotate_BAM(inbam_ref,inbam_annot,outbam, exons_gff, outbam_noannot=None, annot_file=None, strand=1, ambivalence_db=None, ambivalence_db_GENE=None, gene_level_ambivalence=True, outfile_stats=None, updateflags=True, maxd=200, maxd_genomic=21, nmax=0, ambivalence_file_out=None, keep_noannot=True, collapse_annotations=True):
    annot_dict=None
    
    # Load annotation dictionnary
    if not(annot_file==None):
        if isinstance(annot_file,str):
            print("Loading annotations dictionnary from %s"%annot_file)
            annot_dict=gtfutils.make_tx_dict_fromCSV(annot_file)
        else:
            print('Annotation dictionnary passed directly as argument')
            annot_dict=annot_file
        
    # Load exons dictionnary
    exons_dict=None
    if isinstance(exons_gff,str):
        if exons_gff.endswith(".pickle"):
            print("Loading exons dictionnary from pickle file %s "%exons_gff)
            with open(exons_gff, 'rb') as handle:
                exons_dict=pickle.load(handle)
            # exons_dict = gtfutils.load_dict_fromPickle(exons_gff)
        else:
            print("Creating exons dictionnary from %s "%exons_gff)
            exons_dict=gtfutils.make_exons_dict_fromGFF(exons_gff)
    else:
        exons_dict=exons_gff

    # Prepare the ambilvalence dictionnary unless one was passed 
    if not ambivalence_db is None:
        if isinstance(ambivalence_db, str):
            print("Loading ambivalence dictionnary from %s"%ambivalence_db)
            with open(ambivalence_db, 'wb') as handle:
                ambivalence_dict=pickle.load(handle)
        else:
            print('Ambivalence dictionnary passed directly as argument')
            ambivalence_dict=ambivalence_db
        ixmax=max(v for _, v in ambivalence_db.items())
        print('HEY')
        print(ixmax)
    else:
        print("Creating new ambivalence dictionnary")
        ambivalence_dict={'ENST*':1}
        ixmax=1

    
    if not ambivalence_db_GENE is None:
        if isinstance(ambivalence_db_GENE, str):
            print("Loading gene ambivalence dictionnary from %s"%ambivalence_db_GENE)
            with open(ambivalence_db, 'wb') as handle:
                ambivalence_dict_GENE=pickle.load(handle)
        else:
            print('Gene ambivalence dictionnary passed directly as argument')
            ambivalence_dict_GENE=ambivalence_db_GENE
        ixmax_GENE=max(v for _, v in ambivalence_db_GENE.items())
        print(ixmax_GENE)
    else:
        if gene_level_ambivalence:
            print("Creating new gene ambivalence dictionnary")
            ambivalence_dict_GENE={'ENSG*':1}
            ixmax_GENE=1
        else:
            ambivalence_dict_GENE=None
            ixmax_GENE=1
    
    print("Started processing bams")
    isunstranded=False
    isrevstrand=False
    if strand==0:
        isunstranded=True
    if strand==-1:
        isrevstrand=True

    in_ref = pysam.AlignmentFile(inbam_ref)
    in_annot = pysam.AlignmentFile(inbam_annot)
    bamiter_ref = in_ref.fetch(until_eof=True)
    bamiter_annot = in_annot.fetch(until_eof=True)
    outmode="wb" if outbam.endswith(".bam") else "w"
    bam_out = pysam.AlignmentFile(outbam,outmode, template=in_ref) #add header
    
    if keep_noannot and (not outbam_noannot is None):
        bamout_noannot = pysam.AlignmentFile(outbam_noannot,outmode, template=in_ref)
    else:
        bamout_noannot=bam_out

    d_tot, d_stats=_annotate_BAM_stream(bamiter_ref,bamiter_annot,exons_dict, annot_dict, ambivalence_dict, ixmax, ambivalence_dict_GENE, ixmax_GENE, bam_out, bamout_noannot, 
                                isunstranded, isrevstrand, updateflags, keep_noannot, collapse_annotations=collapse_annotations, maxd=maxd, maxd_genomic=maxd_genomic, nmax=nmax)

    if not(outfile_stats==None):
        OUT_stats=open(outfile_stats, "w")
        scalars2report=[k for k, _ in d_tot.items()]
       
        OUT_stats.write(">>TOTALS\n")
        OUT_stats.write(",".join(scalars2report))
        OUT_stats.write("\n")
        OUT_stats.write(",".join([str(d_tot[s2r]) for s2r in scalars2report]))

        dist2report_genomic=['stats_gmappings','stats_gmappings_annotSENSE','stats_gmappings_annotANTISENSE','stats_gmappings_annotNONE',"stats_gmappings_unanottatedONLY"]
        OUT_stats.write("\n\n>>COUNTS_GENOMIC\n")
        OUT_stats.write(",".join(["N"]+dist2report_genomic))
        OUT_stats.write("\n")
        for i in range(maxd_genomic+1):
            v=",".join([str(d_stats[d2r][i]) for d2r in dist2report_genomic])
            OUT_stats.write(str(i)+","+v+"\n")

        dist2report_annot=['stats_annot_ANY','stats_annot_SENSE','stats_annot_ANTISENSE']
        OUT_stats.write("\n>>COUNTS_ANNOTATIONS\n")
        OUT_stats.write(",".join(["N"]+dist2report_annot))
        OUT_stats.write("\n")
        for i in range(maxd+1):
            v=",".join([str(d_stats[d2r][i]) for d2r in dist2report_annot])
            OUT_stats.write(str(i)+","+v+"\n")

        OUT_stats.close()

    #close streams
    bam_out.close()
    in_ref.close()
    in_annot.close()
    if keep_noannot and (not outbam_noannot is None):
        bamout_noannot.close()

    if not ambivalence_file_out is None:
        print("Dumping ambivalence classes to json file %s"%ambivalence_file_out)
        with open(ambivalence_file_out, 'w') as fp:
            json.dump([ambivalence_dict,ambivalence_dict_GENE], fp, indent=4)
        
    print('DONE!')
    return d_tot, d_stats, ambivalence_dict, ambivalence_dict_GENE

