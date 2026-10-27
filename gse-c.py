# Run GSE classifier for continuous predictors
# Take 2 arguments : structure and model
# Read data from StdIn
# continuous predictors

import pyagrum.lib.image as gumimage
from GraphSurgery import GraphSurgeryClassifier
from GraphModel import * 
import pandas as pd
import sys
from sklearn.model_selection import train_test_split
import re
import pdb

#ntiles = 4
structure = sys.argv[1]
model = sys.argv[2]
save_image = sys.argv[3]
tau = float(sys.argv[4])
tau_beta = float(sys.argv[5])

if structure == "observed":
    columns = ['x1', 'x2', 'x3', 'x4', 'x5']
elif structure == "unobserved":
    columns = ['x2', 'x3', 'x4', 'x5']
else:
    columns = ['l1', 'x2', 'x3', 'x4', 'x5']
    



### Read input data from Stdin

df = pd.read_csv(sys.stdin)
dpop = df[df['pop'] == 'dpop']
tpop = df[df['pop'] == 'tpop']
#
#dpop = pd.read_csv("~/dpop.csv")
#tpop = pd.read_csv("~/tpop.csv")

X = dpop[columns]
y = dpop['y'].astype(int)
X_target = tpop[columns]
y_target = tpop['y'].astype(int)

#X = X.apply(lambda x: pd.qcut(x, q=ntiles, labels=range(ntiles)))
#X_target = X_target.apply(lambda x: pd.qcut(x, q=ntiles, labels=range(ntiles)))

X_train, X_validate, y_train, y_validate = train_test_split(X, y, train_size = 0.7, random_state=1)

latent = False
if structure == "observed":
    bn_text = f'''x1 -> y; x2 -> y; x3 -> y; x4 -> y; x5 -> y;'''
elif structure == "unobserved":
    bn_text = f'''x2 -> y; x3 -> y; x4 -> y; x5 -> y;'''
elif structure == "proxy":
    bn_text = f'''x2 -> y; x3 -> y; x4 -> y; x5 -> y; l1 <- u1 -> y;'''
    latent = True


bn = gum.fastBN(bn_text)
target = 'y'

## Specify hidden variable
if latent:
    ## latent_spec (list): a list of tuple that contains the specification of the each latent variable to the observed variables
    latent_spec=[('u1', ('y', 'l1' ))]
    ## a dictionary that maps each variable name to all of its latent confounder as a string to list mapping
    latent_conf = [{'l1': ['u1'], 'y': ['u1']}]
else:
    latent_spec = []
    latent_conf = []

# children of confounding
conf_children = {}
for i in latent_spec:
    latent_name, tup = i
    conf_children[latent_name] = [tup[0], tup[1]]

# Create BN Graph
g = GraphModel()
g.bn_with_latent_as_observed_without_S = bn
g.latent_spec = latent_spec
g.conf_child = conf_children
g.latent_conf = latent_conf 
g.bn_without_latent = g.removeLatent(bn, conf_children)
g.target = target

# construct the bn without u by changing every bidirected edge to directed edge 
g.create_bn_without_latent()

if tau > 0:
    ######################
    # covariate shift
    ######################
    shift = 'cov'
    # Add a selection variable
    if structure == "observed":
        if model == "GSE1":
            chS = ["x1"]
        elif model == "GSE2":
            chS = ["x2"]
    elif structure == "unobserved":
        if model == "GSE1":
            chS = ["y"]
        if model == "GSE2":
            chS = ["x2"]
    elif structure == "proxy":
        if model == "GSE1":
            chS = ['l1', 'y']
        elif model == "GSE2":
            chS = ['x2']

elif tau_beta > 0:
    shift = 'pd'
    ######################
    # pd shift
    ######################
    # Add a selection variable
    if structure == "observed":
        if model == "GSE1":
            chS = ["x1", "y"]
        elif model == "GSE2":
            chS = ["x2"]
    elif structure == "unobserved":
        if model == "GSE1":
            chS = ["y"]
        if model == "GSE2":
            chS = ["x2"]
    elif structure == "proxy":
        if model == "GSE1":
            chS = ['y']
        elif model == "GSE2":
            chS = ['x2']

elif tau_beta == 0 and tau == 0:
    shift = "no"
    ######################
    # pd shift
    ######################
    # Add a selection variable
    if structure == "observed":
        if model == "GSE1":
            chS = []
        elif model == "GSE2":
            chS = ["x2"]
    elif structure == "unobserved":
        if model == "GSE1":
            chS = []
        if model == "GSE2":
            chS = ["x2"]
    elif structure == "proxy":
        if model == "GSE1":
            chS = []
        elif model == "GSE2":
            chS = ['x2']

g.chS = chS
g.s_name_list = ["S"]
bn_with_S = gum.BayesNet(g.bn_without_latent)
bn_with_S.add("S", 2)

# Loop through all nodes in chS and add an arc from S to each one
for node in chS:
    bn_with_S.addArc("S", node)

g.bn_with_S_and_without_latent = bn_with_S
if g.added_edges_tuple:
    g.remove_added_directed_edges_from_bn_without_latent()

    
#remove y -> l1 arc
if structure == "proxy":
    if bn_with_S.existsArc("l1", "y"):
        bn_with_S.eraseArc("l1", "y")

    # convert bn to causal models
g.convertBayesNetToCausalModel()

# Save image
if save_image == 'T':
    gumimage.export(g.cm_with_S, f'~/simulation/gse/{shift}_{structure}_{model}.png')

#Graph surgery 
time_restriction = 120
gse = GraphSurgeryClassifier(g)
gse.fit(X_train, y_train, X_validate, y_validate, zero_one_loss, time_restriction)

if gse.failed:
    print( pd.DataFrame({"predicted_gse": [pd.NA] * tpop.shape[0],"expr": ["failed"] * tpop.shape[0]}).to_csv(index=False) )
else:
    conditioning_set = gse.conditioning_set
    intervention_set = gse.intervention_set
    conditioning_set.sort()
    intervention_set.sort()
    expr = str(gse.id_result)
    full_expr = f'P({target}|do({",".join(intervention_set)}),{",".join(conditioning_set)}) = {expr}'

    # predict target
    y_target_pred = gse.predict_proba_mc(gse.id_result, X, y, X_target, y_target).to_frame()
    y_target_pred['expr'] = full_expr
    #y_target_pred['intervention_conditioning_set'] = intervention_conditioning_set
    #pdb.set_trace()
    #results = list(f'list(y_target_pred = "{y_target_pred}", expr = "{expr}"), intervention_conditioning_set = "{intervention_conditioning_set}"')

    print(y_target_pred.to_csv(index=False))

#print(y_target_pred[["y_target_pred", "expr", "intervention_conditioning_set"]].to_csv(index=False))
