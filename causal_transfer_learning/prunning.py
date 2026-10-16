import pandas as pd
import sys
sys.path.append('causal_transfer_learning/code')
import subset_search

### Read input data

#df = pd.read_csv(sys.argv[1])

### Read input data from Stdin
df = pd.read_csv(sys.stdin)
#df = pd.read_csv("~/dpop.csv")
### Read X
### Format = x1+x2+x3...
X_cols = sys.argv[1].split("+")

### Read y
y_col = sys.argv[2]

X = df[X_cols]
y = df[[y_col]]


clusters = df['cluster'].drop_duplicates().shape[0]
N = df.groupby('cluster').size().iloc[0]
n_ex = [N for i in range(clusters)]
#n_ex = [X.shape[0]]


s_hat = subset_search.subset(X, y, n_ex, valid_split=0.7, delta=0.05, use_hsic=True)

# use greedy search instead to reduce computational time
### s_hat = subset_search.greedy_subset(X, y, n_ex, valid_split=0.8, delta=0.05, use_hsic=True)


# if s_hat returns none of predictors, set formular to 1

if len(s_hat) == 0:
    print(y_col + " ~ 1")
else:
    ### return formul
    print(y_col + " ~ " + " + ".join(X.columns[s_hat]))


### return subset
#print("c" + str(tuple(s_hat+1)))
