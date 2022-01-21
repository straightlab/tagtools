from ..quant import Quantifier

def add_subcommand_quant(subparsers):
    parser = subparsers.add_parser('quant')
    parser.set_defaults(func=quant)
    parser.add_argument('--version', action='version', version='1.0.0')
    parser.add_argument('inbam', help='tagged bam to quantifiy. Should be in compressed form')
    parser.add_argument('tagfile', help='annotation bam')
    parser.add_argument('expressionfile', help='exons description')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--simplify', help='simplify output', action='store_true')
    parser.add_argument('--drop_zeros', help='drop zeros', action='store_true')


    # parser.set_defaults(func=tag)

def quant(args):
    print('resolve', args)
    q=Quantifier(args.tagfile, args.expressionfile)

    
    df=q.quant(args.inbam, todf=True, simplify=args.simplify, drop_zeros=args.drop_zeros)
    df.to_csv(args.outfile)
    