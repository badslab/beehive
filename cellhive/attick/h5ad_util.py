

import pandas as pd
    

def coc(adata, col):
    print(col)
    c = adata.obs[col]
    print('dtype  | ', c.dtype)
    print('head   | ', ' '.join(map(str, c.head())))

def toint(adata, col):
    print(col)
    c = adata.obs[col]
    print('dtype  | ', c.dtype)
    print('head   | ', ' '.join(map(str, c.head())))
    cc = c.astype(int)
    print('head   | ', ' '.join(map(str, cc.head())))
    adata.obs[col] = cc

def tocat(adata, col):
    print(col)
    c = adata.obs[col]
    print('dtype  | ', c.dtype)
    print('head   | ', ' '.join(map(str, c.head())))
    cc = c.astype('category')
    print('head   | ', ' '.join(map(str, cc.head())))
    adata.obs[col] = cc

def toicat(adata, col):
    print(col)
    c = adata.obs[col]
    print('dtype  | ', c.dtype)
    print('head   | ', ' '.join(map(str, c.head())))
    cc = c.astype('int').astype('category')
    print('head   | ', ' '.join(map(str, cc.head())))
    adata.obs[col] = cc


def tofloat(adata, col, onfail=None):
    print(col)
    c = adata.obs[col]
    print('dtype  | ', c.dtype)
    print('head   | ', ' '.join(map(str, c.head())))
    try:
        cc = c.astype('int').astype('float')
    except:
        failed = []
        def force(v, onfail=onfail):
            try:
                return float(v)
            except:
                failed.append(v)
                return onfail

        cc = c.apply(force).astype('float')
        
        if onfail is None:    
            print("Can't convert:")
            print(pd.Series(failed).value_counts())
            return

    print('head   | ', ' '.join(map(str, cc.head())))
    adata.obs[col] = cc

    
def overview(adata):
    import pandas as pd
    r = []
    for n in adata.obs:
        c = adata.obs[n]
        r.append(dict(
            name=n,
            dtype=str(c.dtype),
            head=' | '.join(map(str, c.head()))))
    rv = pd.DataFrame(r)
    print(rv)
    
