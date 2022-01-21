def make_ambiv_LUT(ambivalence_db):
     #db as "ENSTxxx|ENSTxx: 2, etc..."
    #returns "ENSTxx: [2,10,18], ENSTxx: [5], etc"
    # if isinstance(ambivalence_db, str):
    #     print("Loading ambivalence dictionnary from %s"%ambivalence_db)
    #     with open(ambivalence_db, 'wb') as handle:
    #         ambivalence_dict=pickle.load(handle)
    # else:
    #     print('Ambivalence dictionnary passed directly as argument')
    #     ambivalence_dict=ambivalence_db
    
    n=max((v for _, v in ambivalence_db.items()))
    ambiv_LUT={}
    for k, v in ambivalence_db.items():
        tx_all=k.split("|")
        for tx in tx_all:
            tx_ambiv=ambiv_LUT.get(tx, None)
            if tx_ambiv is None:
                ambiv_LUT[tx]=[int(v)]
            else:
                ambiv_LUT[tx]=(tx_ambiv+[int(v)])
    return ambiv_LUT, n

def get_ambiv_truthtable(ambivalence_LUT, nmax, tx_list):
    # [true, false, false, true, ...] ixmax long
    isingroup_truthtable=[False]*nmax
    for tx in tx_list:
        gs=ambivalence_LUT.get(tx, [0])
        for g in gs:
            isingroup_truthtable[g]|=True
    return isingroup_truthtable

