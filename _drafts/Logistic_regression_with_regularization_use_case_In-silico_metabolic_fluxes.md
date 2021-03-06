
---
layout: post
title:  "Logistic Regression with regularization"
date:   2019-03-05
categories: Statistics
---

In one of my previous blogs, I mentioned a technique "FBA", an in-silico method to generate metabolic fluxes. The method has multiple advantages to quickly and efficiently simulate how a model organism's phenotype looks like, given certain input parameters such as the in-silico media on which the organism is fed. Knocking out genes, changing the concentration of the metabolites in media, introducing new reations to the model are other advantages. Refer to the original FBA paper <cite> for more information.


However, here I focus on a previously published work from Wilke lab (I was part of it as well!). More info on the paper can be found here: https://www.ncbi.nlm.nih.gov/pubmed/25502413. Instead of going through the entire work already published, I will talk about a section where the goal is to identify the growth conditions, given the metabolic in-silico flux. So, this is kind of a inverse problem i.e., instead of generating fluxes using FBA, we use fluxes to *predict* growth conditions using logistic regression (with growth conditions as labels).

In the original publication, we used R for statistical analyses (GLMNET package). Here, I use logistic regression module from scikit-learn.

Here is the link for the R script:
https://github.com/clauswilke/Ecoli_FBA_input_prediction/blob/master/Analysis/Scripts/GLMNET.R

The goal now is to write a script that does *almost* similar thing, but using python.

```python
from __future__ import print_function
import cobra
# import deep learning and other modules
import numpy as np
import pandas as pd
from os.path import join

```


```python
# fix random seed for reproducibility
seed = 7
np.random.seed(7)
```


```python
# Flux dataset from predicting bacterial growth conditions study -- for current purposes, THIS IS MOSTLY CONSIDERED RANDOM SYNTHETIC DATA
dataset1 = pd.read_csv("../Data/syntheticFluxData.csv", delimiter=',', header=None)
```


```python
dataset1.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>7</th>
      <th>8</th>
      <th>9</th>
      <th>...</th>
      <th>2375</th>
      <th>2376</th>
      <th>2377</th>
      <th>2378</th>
      <th>2379</th>
      <th>2380</th>
      <th>2381</th>
      <th>2382</th>
      <th>2383</th>
      <th>2384</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>-0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0.004649</td>
      <td>0</td>
      <td>0.004649</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.000000e+00</td>
      <td>4.692800e-34</td>
      <td>1.025800e-27</td>
      <td>6.151300e-29</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>-0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0.004571</td>
      <td>0</td>
      <td>0.004571</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1.348900e-29</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>7.639300e-28</td>
      <td>1.029200e-27</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>-0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0.006978</td>
      <td>0</td>
      <td>0.006978</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>-0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0.004584</td>
      <td>0</td>
      <td>0.004584</td>
      <td>1</td>
      <td>4</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>2.695900e-29</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>0.0</td>
      <td>-0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0.004760</td>
      <td>0</td>
      <td>0.004760</td>
      <td>1</td>
      <td>5</td>
      <td>5</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 2385 columns</p>
</div>




```python
# First 2381 columns are the in-silico fluxes generated from the flux balance analyses (FBA) for the input growth conditions
# Information on the growth conditions is in columns 2382 and 2383. Column 2384 is just the pair-wise combination of 2382 and 2383.
# Original work has lot of details, and this is just a tutorial on how the processed FBA data looks and using deep learning models to analyze.
# Focus is on deep learning models, rather than the source or content of the data.
print("Background information of the data and what we are trying to accomplish here - after importing modules and loading data :-)")

```

    Background information of the data and what we are trying to accomplish here - after importing modules and loading data :-)



```python
# Assuming iAF1260 model used in this repo has exchange/transport reactions that match to the synthetic data.
sbml_path = join("../Data","iAF1260.xml.gz")
iAF1260_ecoli_model = cobra.io.read_sbml_model(sbml_path)
iAF1260_reaction_IDs = [x.id for x in iAF1260_ecoli_model.reactions]
# There might be a better way to do this, but this explicitly conveys information
list_tpp_tex = []
for x in iAF1260_reaction_IDs:
    if(x.endswith('tpp')):
        list_tpp_tex.append(x)
    if(x.endswith('tex')):
        list_tpp_tex.append(x)
# columns to remove from above dataset
IDs_of_cols_to_remove = [iAF1260_reaction_IDs.index(i) for i in list_tpp_tex]
# append output columns as well that are to be removed from Input
IDs_of_cols_to_remove.extend(range(2382,2385)) # Note range covers 2382 to 2384

```


```python
# retain columns that are NOT IDs_of_cols_to_remove, fit in regression methods automatically scales these
X = dataset1.iloc[:,dataset1.columns.difference(IDs_of_cols_to_remove)]
# Use carbon source as y (i.e., col 2382) --- Nitrogen source is 2383, while pair-wise C/N is 2384
fba_y = dataset1.iloc[:,2382]
```


```python
# Two ways to approach this multi-class problem using scikit-learn Logistic regrssion models are 1. ovr, 2. multinomial
# The goal is to predict the growth condition given the metabolic fluxes, so I think may be "ovr" i.e., one vs rest model is better.
# Let's try that for now :-) - If we have to go back and re-do with multinomial and compare the results, we can always do that.
print("Background information of the data")
```

    Background information of the data



```python
X.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
      <th>11</th>
      <th>12</th>
      <th>17</th>
      <th>...</th>
      <th>2368</th>
      <th>2369</th>
      <th>2370</th>
      <th>2371</th>
      <th>2372</th>
      <th>2374</th>
      <th>2375</th>
      <th>2377</th>
      <th>2378</th>
      <th>2380</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.2</td>
      <td>0.0</td>
      <td>...</td>
      <td>-0.0</td>
      <td>0.6</td>
      <td>-0.0</td>
      <td>-0.0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.000000e+00</td>
      <td>4.692800e-34</td>
      <td>1.025800e-27</td>
      <td>6.151300e-29</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>-0.0</td>
      <td>0.8</td>
      <td>-0.0</td>
      <td>0.2</td>
      <td>0.2</td>
      <td>0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1.348900e-29</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>7.639300e-28</td>
      <td>1.029200e-27</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.2</td>
      <td>0.0</td>
      <td>...</td>
      <td>-0.0</td>
      <td>0.6</td>
      <td>-0.0</td>
      <td>-0.0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>-0.0</td>
      <td>0.6</td>
      <td>-0.0</td>
      <td>-0.0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>0.000000e+00</td>
      <td>2.695900e-29</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>...</td>
      <td>-0.0</td>
      <td>0.6</td>
      <td>-0.0</td>
      <td>-0.0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0.0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 2068 columns</p>
</div>




```python
# Split data into training and test (note the encoded_y is a single column)
from sklearn.model_selection import train_test_split
X_train, X_test, fba_y_train, fba_y_test = train_test_split(X, fba_y, test_size=0.25, random_state=0)
```


```python
from sklearn.linear_model import Lasso
lasso_fba = Lasso(alpha=1)
```


```python
lasso_fba.fit(X_train, fba_y_train)
```




    Lasso(alpha=1, copy_X=True, fit_intercept=True, max_iter=1000,
       normalize=False, positive=False, precompute=False, random_state=None,
       selection='cyclic', tol=0.0001, warm_start=False)




```python
# Let's see out of couple of thousand features, which are at least greater than 0.01
lasso_fba.coef_[lasso_fba.coef_ > 0.01]
```




    array([0.13498902, 0.05082819, 0.01110022, 0.08149858])




```python
# Use test data to predict
y_predict = lasso_fba.predict(X_test)
```


```python
## Looks like we did liinear regresssion, not  logistic regression ( a simple test !)
import matplotlib.pyplot as plt
plt.scatter(y_predict, fba_y_test)
plt.xlabel('y-predict')
plt.ylabel('y-test')

```




    Text(0, 0.5, 'y-test')




![png](output_14_1.png)



```python
# THIS IS THE KEY FOR THE ENTIRE SCRIPT - LOOK AT THE PENALTY WHICH IS POSSIBLE WITH SOLVER AND ONE-VERSUS-REST METHOD
# ALSO NOTE THE NUMBER OF ITERATIONS THAT ARE VERY VERY HIGH -- I USED THE SAME ITERATIONS FROM THE ORIGINAL PAPER
# LOWER ITERATIONS MIGHT NOT GET CONVERGED!
from sklearn.linear_model import LogisticRegression
log_lasso_fba = LogisticRegression(penalty='l2', solver='sag', multi_class='ovr', max_iter=9000000)
```


```python
log_lasso_fba.fit(X_train, fba_y_train)
```




    LogisticRegression(C=1.0, class_weight=None, dual=False, fit_intercept=True,
              intercept_scaling=1, max_iter=9000000, multi_class='ovr',
              n_jobs=None, penalty='l2', random_state=None, solver='sag',
              tol=0.0001, verbose=0, warm_start=False)




```python
fba_y_test_log = log_lasso_fba.predict(X_test)
```




#Original R script can be found here: https://github.com/clauswilke/Ecoli_FBA_input_prediction/blob/master/Analysis/Scripts/GLMNET.R (from Wilke lab)



```python
# confusion matrix
from sklearn.metrics import confusion_matrix
confusion_matrix(fba_y_test, fba_y_test_log)
```




    array([[143,   0,   0,  25,   0,   1,   0],
           [  0, 151,   0,  22,   0,   2,   0],
           [  0,   0, 157,  18,   0,   0,   0],
           [  0,   0,   0, 177,   1,   0,   0],
           [  0,   0,   0,  32, 138,   1,   0],
           [  0,   0,   0,  29,   0, 146,   0],
           [  0,   0,   2,  29,   0,   0, 151]], dtype=int64)





I have to go back and see what 4th column is. This growth condition gets predicted as other types as well.
https://github.com/clauswilke/Ecoli_FBA_input_prediction/blob/master/Manuscript/Figures/Fig3.pdf
Same as the case from the original study where other sources sometimes get predicted as acetate - for detailed, explanation read the original paper from Wilke lab.


```python
data_4_heatmap = confusion_matrix(fba_y_test, fba_y_test_log)
import matplotlib.pyplot as plt
import seaborn as sns
ax = sns.heatmap(data_4_heatmap, linewidth=0.5)
plt.xlabel('y-predict')
plt.ylabel('y-test')
plt.show()
```


![png](output_21_0.png)
