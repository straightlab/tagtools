from ..tagger import annotate_BAM
import json

def add_subcommand_tag(subparsers):
    parser = subparsers.add_parser('tag')
    parser.set_defaults(func=tag)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('inbam_ref', help='reference bam')
    parser.add_argument('inbam_annot', help='annotation bam')
    parser.add_argument('exons_gff', help='exons description')
    parser.add_argument('output', help='output file')
    parser.add_argument('--outnoannot', help='output bam for unanottated alignments')
    parser.add_argument('--annotfile', help='gene annotations')
    parser.add_argument('--outstats', help='output statistics file')
    parser.add_argument('--ambivout', help='ambivalence output file')
    parser.add_argument('--addtoambiv', help='une existing ambivalence file', action='store_true') 
    parser.add_argument('--strand', help='strand', type=int, choices=[-1,0,1])


    # parser.set_defaults(func=tag)

def tag(args):
    print('tag', args)
    if args.addtoambiv:
        try:
            with open(args.ambivout, "r") as read_file:
                ambivalence_db_all = json.load(read_file)
                ambivalence_db=ambivalence_db_all[0]
                ambivalence_db_GENE=ambivalence_db_all[1]
        except FileNotFoundError:
            ambivalence_db=None
            ambivalence_db_GENE=None
    else:
        ambivalence_db=None
        ambivalence_db_GENE=None
    
    annotate_BAM(args.inbam_ref,args.inbam_annot,args.output, args.exons_gff, outbam_noannot=args.outnoannot, annot_file=args.annotfile, strand=args.strand, ambivalence_db=ambivalence_db, ambivalence_db_GENE=ambivalence_db_GENE, gene_level_ambivalence=True, outfile_stats=args.outstats, updateflags=True, maxd=200, maxd_genomic=21, nmax=0, ambivalence_file_out=args.ambivout, keep_noannot=True, collapse_annotations=True)

