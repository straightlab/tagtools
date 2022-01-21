#!/usr/bin/env python
# import os
import sys
import pysam
import itertools

def bam2pairs_record(rR,rD,nr,nd,invertFlag,triangular_mode):

    is_reverseD=rD.is_reverse^invertFlag
    is_reverseR=rR.is_reverse^invertFlag
    # s+=" "+"-"*x+ "+"*(not(x))
    # mapq, gap2bridge, AS, XS, secondary_flag
    if is_reverseR:
        posR=rR.reference_start #5'end of RNA
        sR="-"

    else:
        posR=rR.reference_end-1 #5'end of DNA
        sR="+"
        # lenR=str(-rR.reference_length)
    if is_reverseD:
        posD=rD.reference_end-1 #5'end of RNA
        sD="-"
    else:
        posD=rD.reference_start #5'end of DNA
        sD="+"

    gap2bridgeD= str((rD.query_length-rD.query_alignment_end) if invertFlag else rD.query_alignment_start)
    gap2bridgeR= str(rR.query_alignment_start if invertFlag else (rR.query_length-rR.query_alignment_end))

    # xsR=str(rR.get_tag('XS') if rR.has_tag('XS') else 0)
    # xsD=str(rD.get_tag('XS') if rD.has_tag('XS') else 0)
    # asR=str(rR.get_tag('AS') if rR.has_tag('AS') else 0)
    # asD=str(rD.get_tag('AS') if rD.has_tag('AS') else 0)
    
    annot_ref=(rR.get_tag('ar') if rR.has_tag('ar') else "*")
    annot_pos=str(rR.get_tag('ap') if rR.has_tag('ap') else 0)
    annot_name=(rR.get_tag('an') if rR.has_tag('an') else "*")
    annot_type=(rR.get_tag('at') if rR.has_tag('at') else "*")
    annot_strand=(rR.get_tag('as') if rR.has_tag('as') else "*")
    annot_length=(str(rR.get_tag('al')) if rR.has_tag('al') else "0")

    annot_nU=str(rR.get_tag('aU') if rR.has_tag('aU') else 0)
    annot_nM=str(rR.get_tag('aM') if rR.has_tag('aM') else 0)
    annot_nP=str(rR.get_tag('aP') if rR.has_tag('aP') else 0)
    
    nR=str(rR.get_tag('NH') if rR.has_tag('NH') else nr)
    nD=str(rD.get_tag('NH') if rD.has_tag('NH') else nd)

    # is_lowtrig=_lt_chrpos(rD.reference_name,rR.reference_name,int(posD),int(posR))
    # TODO : CATCH bad reference ids when unknown target
    is_lowtrig=_lt_chrpos(rD.reference_id,rR.reference_id,int(posD),int(posR))
    # refD=rD.reference_name
    refR=(triangular_mode==1)*"R_"+rR.reference_name #add a prefix to the RNA which makes it always be smaller than the DNA
    s=""
    if ((triangular_mode==2) & is_lowtrig): #reverse the R and D, do not add prefix
        # readID, chrR, posR, chrD, posD,
        # print("\t".join([rD.query_name,rD.reference_name,str(posD),rR.reference_name,str(posR)]))
        # s="\t".join([rD.query_name,rD.reference_name,str(posD),rR.reference_name,str(posR),sD,sR,str(rD.mapping_quality),str(rR.mapping_quality),asD,asR,xsD,xsR,gap2bridgeD,gap2bridgeR,str(int(rD.is_secondary)),str(int(rR.is_secondary)),str(rD.reference_start),str(rD.reference_end),str(rR.reference_start),str(rR.reference_end),str(nD),str(nR),str(i),annot_ref,annot_pos,annot_name,annot_type, annot_strand, annot_length, str(rD.query_alignment_length), str(rR.query_alignment_length)])
        s="\t".join([
            rD.query_name,rD.reference_name,str(posD),rR.reference_name,str(posR),sD,sR,
            str(rR.mapping_quality),str(rD.mapping_quality), str(rR.flag), str(rD.flag), gap2bridgeR,gap2bridgeD, str(rR.query_alignment_length), str(rD.query_alignment_length),
            nR, nD, annot_nP, annot_nM, annot_nU, annot_ref, annot_pos, annot_name, annot_type, annot_strand, annot_length])
    else:
        #OLD OUTPUT 0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 7=QR, 8=QD, 9=asR, 10=asD, 11=xsR, 12=xsD, 13=gapR, 14=gapD, 15=isSecondaryR, 16=isSecondaryD, 17=rna_start, 18=rna_end, 19=dna_start, 20=dna_end, 21=nR, 22=nD, 23=i, 24=annot_ref, 25=annot_name, 26=annot_type, 27=annot_strand, 28=annot_length, 29=query alignment length RNA, 30=query aligment length DNA
        # s="\t".join([rD.query_name,refR,str(posR),rD.reference_name,str(posD),sR,sD,str(rR.mapping_quality),str(rD.mapping_quality),asR,asD,xsR,xsD,gap2bridgeR,gap2bridgeD,str(int(rR.is_secondary)),str(int(rD.is_secondary)),str(rR.reference_start),str(rR.reference_end),str(rD.reference_start),str(rD.reference_end), str(nR), str(nD), str(i), annot_ref, annot_pos, annot_name, annot_type, annot_strand, annot_length, str(rR.query_alignment_length), str(rD.query_alignment_length)])
        
        #NEW OUTPUT 0=readID, 1=chrR, 2=posR, 3=chrD, 4=posD, 5=strandR, 6=strandD, 
        # 7=QR, 8=QD, 9=flagR, 10=flagD, 11=gapR, 12=gapD, 13=query alignment length RNA, 14=query aligment length DNA, 
        # 15=nR, 16=nD, 17=nP, 18=nM, 19=nU 20=annot_ref, 21=annot_pos, 22=annot_name, 23=annot_type, 24=annot_strand, 25=annot_length, 
        s="\t".join([
            rD.query_name,refR,str(posR),rD.reference_name,str(posD),sR,sD,
            str(rR.mapping_quality),str(rD.mapping_quality), str(rR.flag), str(rD.flag), gap2bridgeR,gap2bridgeD, str(rR.query_alignment_length), str(rD.query_alignment_length),
            nR, nD, annot_nP, annot_nM, annot_nU, annot_ref, annot_pos, annot_name, annot_type, annot_strand, annot_length])

       
    return s, is_lowtrig

def fetch2nextID_(bam_reader,r): # r is the current record
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

def annotate_read_RNA(rR,invertFlag):
    is_reverseR=rR.is_reverse^invertFlag
    posR=rR.reference_start if is_reverseR else (rR.reference_end-1) #5'end of RNA
    gap2bridgeR= rR.query_alignment_start if invertFlag else (rR.query_length-rR.query_alignment_end)
    alg_len=rR.reference_length
    rR.set_tag('ct','R',value_type='A')
    rR.set_tag('cp',posR,value_type='i')
    rR.set_tag('cg',gap2bridgeR,value_type='i')
    rR.set_tag('ca',alg_len,value_type='i')
    rR.set_tag('co',int(invertFlag),value_type='i')

def annotate_read_DNA(rD,invertFlag):
    is_reverseD=rD.is_reverse^invertFlag
    posD=(rD.reference_end-1) if is_reverseD else rD.reference_start #5'end of RNA
    gap2bridgeD= (rD.query_length-rD.query_alignment_end) if invertFlag else rD.query_alignment_start
    alg_len=rD.reference_length
    rD.set_tag('ct','D',value_type='A')
    rD.set_tag('cp',posD,value_type='i')
    rD.set_tag('cg',gap2bridgeD,value_type='i')
    rD.set_tag('ca',alg_len,value_type='i')
    rD.set_tag('co',int(invertFlag),value_type='i')

# def annotate_read_pair(rR,rD,iR,iD,nR,nD,invertFlag):
#     is_reverseD=rD.is_reverse^invertFlag
#     posD=(rD.reference_end-1) if is_reverseD else rD.reference_start #5'end of RNA
#     gap2bridgeD= (rD.query_length-rD.query_alignment_end) if invertFlag else rD.query_alignment_start
#     alg_len=rD.reference_length
#     rD.set_tag('dp',posD,value_type='i')
#     rD.set_tag('dg',gap2bridgeD,value_type='i')
#     rD.set_tag('da',alg_len,value_type='i')
#     rD.set_tag('do',int(invertFlag),value_type='i')
#     rD.set_tag('dn',int(nD),value_type='i')  
#     rR.set_tag('rp',posD,value_type='i')
#     rR.set_tag('rg',gap2bridgeD,value_type='i')
#     rR.set_tag('ra',alg_len,value_type='i')
#     rR.set_tag('ro',int(invertFlag),value_type='i')
#     rR.set_tag('dn',int(nD),value_type='i')

def bampipe_RNA(bamiterR,out,invertFlag):
    nR=0
    try:
        while True:
            rR = bamiterR.__next__()
            nR+=1
            # annotate_read_RNA(rR,invertFlag)
            out.write(rR)
    except StopIteration: # no RNA, move everything to DNA
        return nR

def bampipe_DNA(bamiterD,out,invertFlag):
    nD=0
    try:
        while True:
            rD = bamiterD.__next__()
            nD+=1
            # annotate_read_DNA(rD,invertFlag)
            out.write(rD)
    except StopIteration: # no RNA, move everything to DNA
        return nD

def _lt_natural(x,y, naturalSort):
    if naturalSort:
        return ([ int(c) if c.isdigit() else c.lower() for c in x.split(':') ] < [ int(c) if c.isdigit() else c.lower() for c in y.split(':') ] )
    else:
        return x<y

# def _lt_chrpos(x,y,xpos,ypos):
#     if _lt_natural(x,y):
#         return True
#     elif x==y:
#         (return True) if x<y else (return False)
#     else:
#         return False

def _lt_chrpos(x,y,xpos,ypos):
    return ( (x<y) | ((x==y) & (xpos<ypos)) )

def bams2pairs_stream(bamiterR,bamiterD,OUT_pair,OUT_r,OUT_d,invert=False, OUT_pairedRNA=None, OUT_pairedDNA=None, triangular_mode=0, OUT_pair_lowtrig=None, naturalSort=False, OUT_fastq_D=None, OUT_fastq_R=None): #0, keep everything intact, 1=enforce upper triangular by prefixng the rna with R_, 2=enforce upper triangular by splitting into to two files (file 2 has RNA and DNA exchanged and therefore stays upper triang )
    # IMPORTANT: the streams need to generate records sorted by readIDs
    npairs=0
    npairs_lowtrig=0
    nD=0
    nR=0
    try:
        rR = bamiterR.__next__()
    except StopIteration: # no RNA, move everything to DNA and done
        nD=bampipe_DNA(bamiterD,OUT_d,invert)
        return (npairs,nR,nD)
    try:
        rD = bamiterD.__next__()
    except StopIteration: # no DNA, move everything to RNA and done
        nR=bampipe_RNA(bamiterR,OUT_r,invert)
        return (npairs,nR,nD)

    while not((rD==None) & (rR==None)):
        if rD==None: # no more DNA data, finish writing everything to RNA file and done
            # annotate_read_RNA(rR, invert) #we don't annotate these reads anymore, while waste time
            OUT_r.write(rR)
            nRpipe=bampipe_RNA(bamiterR,OUT_r,invert)
            nR+=(1+nRpipe)
            return (npairs,nR,nD)
        elif rR==None: # no more RNA data, finish writing everything to DNA file and done
            # annotate_read_DNA(rD, invert)
            OUT_d.write(rD)
            nDpipe=bampipe_DNA(bamiterD,OUT_d,invert)
            nD+=(1+nDpipe)
            return (npairs,nR,nD)
        else:
            while not(rR.query_name==rD.query_name): # no match
                if  _lt_natural(rD.query_name,rR.query_name, naturalSort): #no match and the DNA stream is behind
                    # annotate_read_DNA(rD, invert)
                    OUT_d.write(rD)
                    nD+=1
                    try:
                        rD=bamiterD.__next__()
                    except StopIteration: # no more DNA data, finish writing everything to RNA file and done
                        # annotate_read_RNA(rR,invert)
                        OUT_r.write(rR)
                        nRpipe=bampipe_RNA(bamiterR,OUT_r,invert)
                        nR+=(1+nRpipe)
                        return (npairs,nR,nD)

                else: # no match, and the RNA stream is behind
                    #annotate_read_RNA(rR, invert)
                    OUT_r.write(rR)
                    nR+=1
                    try:
                        rR=bamiterR.__next__()
                    except StopIteration: # no more RNA data, finish writing everything to RNA file and done
                        #annotate_read_DNA(rD,invert)
                        OUT_d.write(rD)
                        nDpipe=bampipe_RNA(bamiterD,OUT_d,invert)
                        nD+=(1+nDpipe)
                        return (npairs,nR,nD)

            # at this point, the R and D streams are pointing to unanalzyed records with matching IDs, so we gather them along with possible duplicates

            recordsD, rD = fetch2nextID_(bamiterD,rD)
            recordsR, rR = fetch2nextID_(bamiterR,rR)
#             n=bam2pairs_record_list(recordsR,recordsD,OUT_pair,invert)
            npairs+=1


            # no more extra bam files
            # if not(OUT_fastq_D==None):
            #     OUT_fastq_D.write("@%s\n" % recordsD[0].query_name)
            #     OUT_fastq_D.write("%s\n+\n" % recordsD[0].get_forward_sequence())
            #     xx=recordsD[0].get_forward_qualities()
            #     for i in range(len(xx)):
            #         xx[i]+=33
            #     OUT_fastq_D.write("%s\n" % xx.tostring().decode('utf-8'))

            #     OUT_fastq_R.write("@%s\n" % recordsR[0].query_name)
            #     OUT_fastq_R.write("%s\n+\n" % recordsR[0].get_forward_sequence())
            #     xx=recordsR[0].get_forward_qualities()
            #     for i in range(len(xx)):
            #         xx[i]+=33
            #     OUT_fastq_R.write("%s\n" % xx.tostring().decode('utf-8'))

            n_recR=len(recordsR)
            n_recD=len(recordsD)
            # if not(OUT_pairedRNA==None):
            #     for recR in recordsR:
            #         annotate_read_RNA(recR, invert) # TO DO: keep track of duplication here for a duplication stats
            #         OUT_pairedRNA.write(recR)

            # if not(OUT_pairedDNA==None):
            #     for recD in recordsD:
            #         annotate_read_DNA(recD, invert)
            #         OUT_pairedDNA.write(recD)

            #TO DO: offer other modes of operation where we don't take combinatorial matches
            if triangular_mode==2:
                for (_, elem) in enumerate(itertools.product(recordsR,recordsD)):
                    s, is_lowtrig = bam2pairs_record(elem[0],elem[1],n_recR,n_recD,invert,triangular_mode)
                    if is_lowtrig:
                        npairs_lowtrig+=1
                        OUT_pair_lowtrig.write(s)
                        OUT_pair_lowtrig.write("\n")
                    else:
                        OUT_pair.write(s)
                        OUT_pair.write("\n")
            else:
                for (_, elem) in enumerate(itertools.product(recordsR,recordsD)):
                    s, is_lowtrig = bam2pairs_record(elem[0],elem[1],n_recR,n_recD,invert,triangular_mode)
                    OUT_pair.write(s)
                    OUT_pair.write("\n")
                    if is_lowtrig:
                        npairs_lowtrig+=1
    return npairs,nR,nD,npairs_lowtrig


def bams2pairs(infile_rna,infile_dna,outfile_pairs,outfile_rna,outfile_dna,triangular_mode=0, outfile_pairs_lowtrig=None, invert=False, outfile_pairedRNA=None, outfile_pairedDNA=None, naturalSort=False, outfastq_DNA=None, outfastq_RNA=None, outfile_stats=None):
    if invert:
        print("Running in reverse bridge mode")
    inr = pysam.AlignmentFile(infile_rna)
    ind = pysam.AlignmentFile(infile_dna)
    bamiterR = inr.fetch(until_eof=True)
    bamiterD = ind.fetch(until_eof=True)
    OUT_pair = open(outfile_pairs, "w", encoding="utf-8")
    outmode_r="wb" if outfile_rna.endswith(".bam") else "w"
    outmode_d="wb" if outfile_dna.endswith(".bam") else "w"
    OUT_r = pysam.AlignmentFile(outfile_rna,outmode_r, template=inr) #add header
    OUT_d = pysam.AlignmentFile(outfile_dna,outmode_d, template=ind)

    OUT_pairedRNA=None
    if not(outfile_pairedRNA==None):
        outmode_pairedRNA="wb" if outfile_pairedRNA.endswith(".bam") else "w"
        OUT_pairedRNA= pysam.AlignmentFile(outfile_pairedRNA, outmode_pairedRNA, template=inr)

    OUT_pairedDNA=None
    if not(outfile_pairedDNA==None):
        outmode_pairedDNA="wb" if outfile_pairedDNA.endswith(".bam") else "w"
        OUT_pairedDNA= pysam.AlignmentFile(outfile_pairedDNA, outmode_pairedDNA, template=ind)

    OUT_pair_lowtrig=None
    if not(outfile_pairs_lowtrig==None):
        OUT_pair_lowtrig= open(outfile_pairs_lowtrig, "w", encoding="utf-8")

    OUT_fastq_D=None
    OUT_fastq_R=None
    if not(outfastq_DNA==None):
        OUT_fastq_D= open(outfastq_DNA, "w")
        OUT_fastq_R= open(outfastq_RNA, "w")

    
    n=bams2pairs_stream(bamiterR,bamiterD,OUT_pair,OUT_r,OUT_d, invert, OUT_pairedRNA, OUT_pairedDNA, triangular_mode, OUT_pair_lowtrig, naturalSort, OUT_fastq_D,OUT_fastq_R)

    if not(outfile_stats==None):
        OUT_stats=open(outfile_stats, "w")
        OUT_stats.write("n,nRNAsOnly,nDNAsOnly,nlowtrig\n")
        OUT_stats.write(",".join([str(n[0]),str(n[1]),str(n[2]),str(n[3])]))
        OUT_stats.write("\n")
        OUT_stats.close()

    #close streams
    OUT_pair.close()
    OUT_r.close()
    OUT_d.close()
    if not(outfile_pairedRNA==None):
        OUT_pairedRNA.close()
    if not(outfile_pairedDNA==None):
        OUT_pairedDNA.close()
    if not(outfile_pairs_lowtrig==None):
        OUT_pair_lowtrig.close()
    if not(outfastq_DNA==None):
        OUT_fastq_D.close()
        OUT_fastq_R.close()
    return n

if __name__ == "__main__":
    pairs_prefix=sys.argv[3][0:-6]
    triangular_mode=int(sys.argv[4])
    naturalSort=bool(sys.argv[5])
    outfile_pairs=pairs_prefix + ".pairs"
    outfile_rna=pairs_prefix + ".unpaired.rna.bam"
    outfile_dna=pairs_prefix + ".unpaired.dna.bam"
    outfile_pairedRNA=None #pairs_prefix + ".paired.rna.bam"
    outfile_pairedDNA=None #pairs_prefix + ".paired.dna.bam"
    outfastq_DNA=None #pairs_prefix + '.dna.fastq'
    outfastq_RNA=None #pairs_prefix + '.rna.fastq'
    outfile_lowtrig_pairs=None
    outfile_stats=pairs_prefix + ".stats.txt"
    if triangular_mode==2:
        outfile_lowtrig_pairs=pairs_prefix + ".lowtrig.pairs"
    bams2pairs(sys.argv[1],sys.argv[2],outfile_pairs,outfile_rna,outfile_dna,triangular_mode=triangular_mode, outfile_pairs_lowtrig=outfile_lowtrig_pairs, invert=False, outfile_pairedRNA=outfile_pairedRNA, outfile_pairedDNA=outfile_pairedDNA, naturalSort=naturalSort,outfastq_DNA=outfastq_DNA, outfastq_RNA=outfastq_RNA, outfile_stats=outfile_stats)
