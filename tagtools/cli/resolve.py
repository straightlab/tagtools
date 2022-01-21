from ..quant import Quantifier

def add_subcommand_resolve(subparsers):
    parser = subparsers.add_parser('resolve')
    parser.set_defaults(func=resolve)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('inbam', help='tagged bam to draw from')
    parser.add_argument('tagfile', help='annotation bam')
    parser.add_argument('expressionfile', help='exons description')
    parser.add_argument('outbam', help='output file')
    parser.add_argument('--annotfile', help='gene annotations')


    # parser.set_defaults(func=tag)

def resolve(args):
    print('resolve', args)
    q=Quantifier(args.tagfile, args.expressionfile)
    if not args.annotfile is None:
        q.set_annotations_dictionary(args.annotfile)
    
    q.resolve(args.inbam, args.outbam)
