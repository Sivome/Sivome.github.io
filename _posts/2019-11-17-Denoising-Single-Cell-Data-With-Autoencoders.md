---
title: "Denoising-Single-Cell-Data-With-Autoencoders"
output: html_document
---



```python
# Load local files to google colab
from google.colab import files
uploaded = files.upload()
```



     <input type="file" id="files-a304b5b9-2904-4c40-bbca-f068d2cc5771" name="files[]" multiple disabled />
     <output id="result-a304b5b9-2904-4c40-bbca-f068d2cc5771">
      Upload widget is only available when the cell has been executed in the
      current browser session. Please rerun this cell to enable.
      </output>
      <script src="/nbextensions/google.colab/files.js"></script>


    Saving mca.10K.vg.csv to mca.10K.vg (1).csv


In my previous blog, I used single cell Mouse Cell Atlas [MCA] data to identify clusters and find differentially expressed markers between the clusters. Single cell datasets, in general are sparse, as the number of genes identified is less, compared to bulk RNA-seq experiment. Another reason for sparsity is the number of cells involved, and in most cases, these cells belong to different cell-types leading to expression of different genes in each cell. There are Python modules / R packages that can denoise the single cell counts data and perform imputation (i.e., identify real zeros to possible artifacts caused by single cell sequencing technologies). Some of the programs that does the denoising + imputation on single cell data are DCA, SAVER, MAGIC, scIMPUTE.

Most of these programs use autoencoders (a variant of neural networks) to learn the structure of the data and fill in the missing values. Autoencoder (as the name suggests) has same number of input and output neurons. In single cell case, if the input is counts data of "m" genes and "n" cells, then the output will have the same dimensions, but has a denoised data with reduced sparsity!

Depending on the number of hidden layers, and activation nodes, there could be different variants of the autoencoders. Here, I wrote a script of an autoencoder implementation that uses additional hidden layers compared to above programs. Another reason I made my own implementation is to understand the basic concepts of how these autoencoders work on single cell data.

The input data I used here is the same I used from my previous blog for Seurat marker analyses i.e., counts data of the publicly available mouse cell atlas. The end goal of this task is to compare the clusters and markers identified with the original counts data, and the autoencoder's denoised data.


```python
import pandas as pd
mca_counts = pd.read_csv("mca.10K.vg.csv", index_col=0)
```


```python
mca_counts.head()
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
      <th>Bladder_1.CGGCAGAAAGTTATTCCA</th>
      <th>Bladder_1.CTCGCAAATAAAATCAAC</th>
      <th>Bladder_1.CCATCTAGCGAGTTTAGG</th>
      <th>Bladder_1.GAGGAGCGCTTGATACAG</th>
      <th>Bladder_1.CCAGACACAATAGAATTA</th>
      <th>Bladder_1.CCGACGGGACATATGGCG</th>
      <th>Bladder_1.TAGCATTCAAAGATTCCA</th>
      <th>Bladder_1.CTCCATCCATCTTTTAGG</th>
      <th>Bladder_1.CAAAGTAGGACTAAGTAC</th>
      <th>Bladder_1.TAGCATTGCGGAAACCTA</th>
      <th>Bladder_1.ACACCCACGTTGCACAAG</th>
      <th>Bladder_1.CGCACCATCTCTTAGTCG</th>
      <th>Bladder_1.ACAATAGTCCCGAACGCC</th>
      <th>Bladder_1.CGAGTATAGCATTCAAAG</th>
      <th>Bladder_1.TCGTAATCGGGTGCAGGA</th>
      <th>Bladder_1.ACAATACCGCTATATGTA</th>
      <th>Bladder_1.CGCACCTGATCATTCATA</th>
      <th>Bladder_1.TTCATAGCGAATGTAATG</th>
      <th>Bladder_1.CGCACCTACTTCTCTACC</th>
      <th>Bladder_1.CAACAAGAGATCAACCTA</th>
      <th>Bladder_1.AGGACTCGCACCGAGATC</th>
      <th>Bladder_1.TAGCATCGTATTGCGAAT</th>
      <th>Bladder_1.CCGACGATGCTTCGCTTG</th>
      <th>Bladder_1.TATGTAAAAACGTCTACC</th>
      <th>Bladder_1.CGTATTGATCTTCCTAGA</th>
      <th>Bladder_1.AACCTACAACAAAAAGTT</th>
      <th>Bladder_1.TCACTTGTCGGTTTCATA</th>
      <th>Bladder_1.TTCATATCAAAGTATTGT</th>
      <th>Bladder_1.ACGAGCGAGATCTATTGT</th>
      <th>Bladder_1.AACCTACCATCTGAGGAG</th>
      <th>Bladder_1.AACGCCAAAGTTCCTTTC</th>
      <th>Bladder_1.ACTTATAGTCGTTTTAGG</th>
      <th>Bladder_1.AAGTACATGGCGGCAGGA</th>
      <th>Bladder_1.CATGATAGTCGTCTTCTG</th>
      <th>Bladder_1.AGCGAGAGCGAGAGATGG</th>
      <th>Bladder_1.ACGTTGGTTGCCCATGAT</th>
      <th>Bladder_1.TGATCACACAAGTCAAAG</th>
      <th>Bladder_1.CCTAGATTCCGCGTATAC</th>
      <th>Bladder_1.TTAACTTCGTAAGCAGGA</th>
      <th>Bladder_1.CCAGACATTCCACTGTGT</th>
      <th>...</th>
      <th>BoneMarrow_5.CCTAGACATGATGCTGTG</th>
      <th>BoneMarrow_5.CGCACCATCTCTTGAAGC</th>
      <th>BoneMarrow_5.TGCAATATCTCTACTTAT</th>
      <th>BoneMarrow_5.ACTTATAAAGTTACAATA</th>
      <th>BoneMarrow_5.CACAAGTGCGGATGCGGA</th>
      <th>BoneMarrow_5.CGTATTATTCCATCTACC</th>
      <th>BoneMarrow_5.TGAAGCCCAGACCTCGCA</th>
      <th>BoneMarrow_5.CCTTTCACTTATCTGTGT</th>
      <th>BoneMarrow_5.CGCTTGCACAAGAGGACT</th>
      <th>BoneMarrow_5.TCACTTAGTTTAGCTGTG</th>
      <th>BoneMarrow_5.ACCTGACGTGGCGTATAC</th>
      <th>BoneMarrow_5.TCAAAGCGGCAGGGGCGA</th>
      <th>BoneMarrow_5.ACACCCCTCCATTCTACC</th>
      <th>BoneMarrow_5.GATCTTGGGCGAATGGCG</th>
      <th>BoneMarrow_5.ACTTATGACACTGGCTGC</th>
      <th>BoneMarrow_5.TTTAGGGGCTGCGGCTGC</th>
      <th>BoneMarrow_5.TGATCATCACTTCCAGAC</th>
      <th>BoneMarrow_5.TGCGGAATTCCAGATCTT</th>
      <th>BoneMarrow_5.CAACAAATCTCTCTGAAA</th>
      <th>BoneMarrow_5.ATCAACTTAACTCCAGAC</th>
      <th>BoneMarrow_5.CTCGCATAGAGAGAACGC</th>
      <th>BoneMarrow_5.ACCTGAAACCTAGCGAAT</th>
      <th>BoneMarrow_5.TGATCACCAGACCGAGTA</th>
      <th>BoneMarrow_5.ACGAGCTACTTCACTTAT</th>
      <th>BoneMarrow_5.GATCTTAGGGTCGGACAT</th>
      <th>BoneMarrow_5.GAGGAGGGACATATGCTT</th>
      <th>BoneMarrow_5.TGATCATGATCAAGGACT</th>
      <th>BoneMarrow_5.AAGCGGAACCTAGCTGTG</th>
      <th>BoneMarrow_5.GAGGAGGGACATTCAAAG</th>
      <th>BoneMarrow_5.TGCAATGATCTTGTCCCG</th>
      <th>BoneMarrow_5.ATGGCGAATAAATGCAAT</th>
      <th>BoneMarrow_5.CCTTTCTTTAGGCTGTGT</th>
      <th>BoneMarrow_5.AGGGTCAACCTATGATCA</th>
      <th>BoneMarrow_5.CTGAAAGCAGGATAGAGA</th>
      <th>BoneMarrow_5.CGCTTGCGAGTAGTCGGT</th>
      <th>BoneMarrow_5.CAACAAGTTGCCTATTGT</th>
      <th>BoneMarrow_5.ACCTGACACAAGTTAACT</th>
      <th>BoneMarrow_5.CGTGGCATTTGCGTAATG</th>
      <th>BoneMarrow_5.GCAGGAGTTGCCGGACAT</th>
      <th>BoneMarrow_5.CTCCATTGTCACCGAGTA</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1110008P14Rik</th>
      <td>6</td>
      <td>2</td>
      <td>6</td>
      <td>5</td>
      <td>7</td>
      <td>8</td>
      <td>7</td>
      <td>6</td>
      <td>5</td>
      <td>8</td>
      <td>7</td>
      <td>9</td>
      <td>7</td>
      <td>10</td>
      <td>6</td>
      <td>3</td>
      <td>6</td>
      <td>5</td>
      <td>5</td>
      <td>8</td>
      <td>3</td>
      <td>5</td>
      <td>5</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>3</td>
      <td>6</td>
      <td>5</td>
      <td>2</td>
      <td>3</td>
      <td>7</td>
      <td>1</td>
      <td>5</td>
      <td>6</td>
      <td>5</td>
      <td>2</td>
      <td>3</td>
      <td>2</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700016P03Rik</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700027J19Rik</th>
      <td>4</td>
      <td>2</td>
      <td>2</td>
      <td>6</td>
      <td>4</td>
      <td>2</td>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>3</td>
      <td>1</td>
      <td>3</td>
      <td>2</td>
      <td>1</td>
      <td>2</td>
      <td>4</td>
      <td>3</td>
      <td>4</td>
      <td>3</td>
      <td>2</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
      <td>7</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700056N10Rik</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700119H24Rik</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 10000 columns</p>
</div>




```python
mca_counts.shape
```




    (2000, 10000)



**bold text**


```python

```

Some information the input dataset ( 2000 variable genes in 10000 single cells ): Instead of using the entire gene lists for which the counts data is available, I used the top 2000 variable genes of the dataset. This implies that the number of training data sets is 2000 (in computer vision applications, this implies the number of images used for training). The shape of the input dimension for images in computer vision applications is x by y (or x by y by 3 in case of RGB). This is synonymous to the number of cells (single-cell) for which gene counts data is available. In a sample MNIST data, i.e., grey scale low-res image, the size of image is 28 by 28 i.e., 784. Here, the size of the data set (also called the input dimension in neural network nomenclature) is equal to the number of cells considered. We used subset of the cells from mouse cell atlas, i.e., 10K cells.

 The input vector of each input neuron is the counts of a particular gene across all cells. In general, counts data are modelled with a poisson distribution (mean and variance are equal) or in most cases, using negative binomial distribution (mean and variance can take different values, compared to poisson). In DCA, and SAVER-X, 2 variants of negative binomial (NB) distribution are implemented, depending on the type of input data used.

The end goal of this blog is to show the basic details of different aspects of denoising single cell data. For accurate classifications and analyses, it is highly recommended to use the above programs, but this blog would serve as reference to understand the basic details of [1] generating the input data for the above tools (like SAVERX, DCA), [2] understanding the notation of the neural networks (i.e., training data set size, input dimension of the input layer neurons), [3] neural network model i.e., layers of the autoencoder i.e., input layer + [ encoder layers ] + [ compression layer ] + [ decoder layer ] + [ output layer ], [4] the parameters used to run the neural network to calculate the loss and other metrics, [5] evaluating the model on the input data to generate the denoised + impute version of the original single cell data. Finally, both the original and the basic denoised version of the dataset is fed back into the Seurat for clustering analyses. PLEASE NOTE that this model lacks using the validation data and the test data (which are mandatory for any ML modeling projects). In my next blog, I can go in depth about using the validation data to develop better ML autoencoder model, and using a different test set.

Other than the model and the optimization functions, loss function should be carefully picked depending on the application. Here, I used poisson loss function that is available in the keras module. The only difference between the previous implementations and the current one is the number of hidden layers used in encoder and decoder. Unlike one layer of encoder and decoder, I implemented 2 layers with 256 and 64 neurons each. The encoded layers have 256 and 64neurons, while the decoded layers has 256 and 64 neurons with a final output layer having 2000 neurons representing the number of features (genes). Except for the final output layer, all the activations used are ReLU. The output layer has an activation of sigmoid to scale the output between [0,1].


```python
from keras.layers import Input, Dense, BatchNormalization, Dropout, Activation
from keras.models import Model
#from keras.regularizers import l1, L1L2, l2

# Input layer
input_cell = Input(shape=(10000, ))

# Encoded layer 1
encoded = Dense(256, activation='relu')(input_cell) #, activity_regularizer=L1L2(l1=10e-3, l2=10e-5)
encoded = Dropout(rate=0.2)(encoded)

# Encoded layer 2
encoded_2 = Dense(64, activation='relu')(encoded)
encoded_2 = Dropout(rate=0.2)(encoded_2)

# Compression layer 2
compression_layer = Dense(16, activation='relu')(encoded_2)
compression_layer = Dropout(rate=0.2)(compression_layer)

# Decoded layer 2
decoded_2 = Dense(64, activation='relu')(compression_layer)
decoded_2 = Dropout(rate=0.2)(decoded_2)

# Decoded layer 1
decoded = Dense(256, activation='relu')(decoded_2)
decoded = Dropout(rate=0.2)(decoded)

# Output layer
output_cell = Dense(10000, activation='sigmoid')(decoded)

# Defining the entire model
autoencoder = Model(input_cell, output_cell)

# compile the model before fit and predict
autoencoder.compile(optimizer='rmsprop', loss='poisson')

# Summary of the model in tabular format to show the  number of trainable and non-trainable parameters
autoencoder.summary()
```

    Model: "model_16"
    _________________________________________________________________
    Layer (type)                 Output Shape              Param #   
    =================================================================
    input_20 (InputLayer)        (None, 10000)             0         
    _________________________________________________________________
    dense_99 (Dense)             (None, 256)               2560256   
    _________________________________________________________________
    dropout_40 (Dropout)         (None, 256)               0         
    _________________________________________________________________
    dense_100 (Dense)            (None, 64)                16448     
    _________________________________________________________________
    dropout_41 (Dropout)         (None, 64)                0         
    _________________________________________________________________
    dense_101 (Dense)            (None, 16)                1040      
    _________________________________________________________________
    dropout_42 (Dropout)         (None, 16)                0         
    _________________________________________________________________
    dense_102 (Dense)            (None, 64)                1088      
    _________________________________________________________________
    dropout_43 (Dropout)         (None, 64)                0         
    _________________________________________________________________
    dense_103 (Dense)            (None, 256)               16640     
    _________________________________________________________________
    dropout_44 (Dropout)         (None, 256)               0         
    _________________________________________________________________
    dense_104 (Dense)            (None, 10000)             2570000   
    =================================================================
    Total params: 5,165,472
    Trainable params: 5,165,472
    Non-trainable params: 0
    _________________________________________________________________



```python
# Log normalizing the raw counts and scaling to 10K cells is typically used with all the single cell programs.
mca_counts_lc = mca_counts.max(axis=0)
mca_counts_norm_lc = mca_counts / mca_counts_lc[None, :]

# Fitting the normalized data to the autoencoder model built
autoencoder.fit(mca_counts_norm_lc, mca_counts_norm_lc, epochs = 20, batch_size = 20, validation_split=0.3)# class_weight=mca_lc_weights, sample_weight=mca_counts_gs)

# Predict the output on the same training data [In real world applications, there should be validation data and a separate test data -- for simplicity, I used the same input data here]
mca_counts_norm_lc_output = autoencoder.predict(mca_counts_norm_lc)

```

    Train on 1400 samples, validate on 600 samples
    Epoch 1/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0187 - val_loss: 0.0294
    Epoch 2/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0188 - val_loss: 0.0290
    Epoch 3/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0187 - val_loss: 0.0294
    Epoch 4/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0189 - val_loss: 0.0298
    Epoch 5/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0185 - val_loss: 0.0289
    Epoch 6/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0186 - val_loss: 0.0298
    Epoch 7/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0187 - val_loss: 0.0301
    Epoch 8/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0187 - val_loss: 0.0302
    Epoch 9/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0189 - val_loss: 0.0298
    Epoch 10/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0188 - val_loss: 0.0297
    Epoch 11/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0185 - val_loss: 0.0299
    Epoch 12/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0188 - val_loss: 0.0296
    Epoch 13/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0185 - val_loss: 0.0300
    Epoch 14/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0184 - val_loss: 0.0291
    Epoch 15/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0188 - val_loss: 0.0291
    Epoch 16/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0186 - val_loss: 0.0295
    Epoch 17/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0185 - val_loss: 0.0294
    Epoch 18/20
    1400/1400 [==============================] - 3s 2ms/step - loss: 0.0186 - val_loss: 0.0293
    Epoch 19/20
    1400/1400 [==============================] - 4s 3ms/step - loss: 0.0185 - val_loss: 0.0293
    Epoch 20/20
    1400/1400 [==============================] - 4s 3ms/step - loss: 0.0185 - val_loss: 0.0297


The loss function used here is the poisson that is available in the keras losses. Another option is the zero inflated poisson model, along with other counts based models such as negative binomial loss functions. The regularization functions (L1, L2 or combination of L1L2 norms), the dropout rate to make sure that the model is not overfit and batch normalization (in case, batch_size is set to greater than 1) are other parameters that can be tuned. The epochs would fine-tune the model, but also take up computation time. So, for the blog purposes, I used the parameters to build the model relatively quickly.

In poisson loss function, the loss would equal to the mean([y_pred - y_true * log(y_pred + epsilon)]), where epsilon is a small number to avoid calculating log(0). Also, for Seurat analyses, it would help to use a threshold and set everything below threshold to zero. However, this threshold should be carefully chosen depending on the dataset (for example, looking at the histogram of the predicted output).


```python
mca_counts_norm_lc.sum(axis=1)
```




    1110008P14Rik    143.055429
    1700016P03Rik      0.124542
    1700027J19Rik     16.903346
    1700056N10Rik      0.155495
    1700119H24Rik      0.222472
                        ...    
    Slc22a3            2.526988
    Treml1             0.255371
    Adamts13           0.089922
    C8a                0.362745
    Ighv5-17           0.131173
    Length: 2000, dtype: float64




```python
mca_counts_norm_lc_output[0:4, 0:4]
```




    array([[0.01773182, 0.02895972, 0.01728064, 0.01625928],
           [0.00141433, 0.00209311, 0.00145644, 0.0011279 ],
           [0.02724421, 0.04175797, 0.02742231, 0.02623108],
           [0.0013988 , 0.00204527, 0.00144061, 0.00110739]], dtype=float32)




```python
mca_counts_norm_lc.head()
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
      <th>Bladder_1.CGGCAGAAAGTTATTCCA</th>
      <th>Bladder_1.CTCGCAAATAAAATCAAC</th>
      <th>Bladder_1.CCATCTAGCGAGTTTAGG</th>
      <th>Bladder_1.GAGGAGCGCTTGATACAG</th>
      <th>Bladder_1.CCAGACACAATAGAATTA</th>
      <th>Bladder_1.CCGACGGGACATATGGCG</th>
      <th>Bladder_1.TAGCATTCAAAGATTCCA</th>
      <th>Bladder_1.CTCCATCCATCTTTTAGG</th>
      <th>Bladder_1.CAAAGTAGGACTAAGTAC</th>
      <th>Bladder_1.TAGCATTGCGGAAACCTA</th>
      <th>Bladder_1.ACACCCACGTTGCACAAG</th>
      <th>Bladder_1.CGCACCATCTCTTAGTCG</th>
      <th>Bladder_1.ACAATAGTCCCGAACGCC</th>
      <th>Bladder_1.CGAGTATAGCATTCAAAG</th>
      <th>Bladder_1.TCGTAATCGGGTGCAGGA</th>
      <th>Bladder_1.ACAATACCGCTATATGTA</th>
      <th>Bladder_1.CGCACCTGATCATTCATA</th>
      <th>Bladder_1.TTCATAGCGAATGTAATG</th>
      <th>Bladder_1.CGCACCTACTTCTCTACC</th>
      <th>Bladder_1.CAACAAGAGATCAACCTA</th>
      <th>Bladder_1.AGGACTCGCACCGAGATC</th>
      <th>Bladder_1.TAGCATCGTATTGCGAAT</th>
      <th>Bladder_1.CCGACGATGCTTCGCTTG</th>
      <th>Bladder_1.TATGTAAAAACGTCTACC</th>
      <th>Bladder_1.CGTATTGATCTTCCTAGA</th>
      <th>Bladder_1.AACCTACAACAAAAAGTT</th>
      <th>Bladder_1.TCACTTGTCGGTTTCATA</th>
      <th>Bladder_1.TTCATATCAAAGTATTGT</th>
      <th>Bladder_1.ACGAGCGAGATCTATTGT</th>
      <th>Bladder_1.AACCTACCATCTGAGGAG</th>
      <th>Bladder_1.AACGCCAAAGTTCCTTTC</th>
      <th>Bladder_1.ACTTATAGTCGTTTTAGG</th>
      <th>Bladder_1.AAGTACATGGCGGCAGGA</th>
      <th>Bladder_1.CATGATAGTCGTCTTCTG</th>
      <th>Bladder_1.AGCGAGAGCGAGAGATGG</th>
      <th>Bladder_1.ACGTTGGTTGCCCATGAT</th>
      <th>Bladder_1.TGATCACACAAGTCAAAG</th>
      <th>Bladder_1.CCTAGATTCCGCGTATAC</th>
      <th>Bladder_1.TTAACTTCGTAAGCAGGA</th>
      <th>Bladder_1.CCAGACATTCCACTGTGT</th>
      <th>...</th>
      <th>BoneMarrow_5.CCTAGACATGATGCTGTG</th>
      <th>BoneMarrow_5.CGCACCATCTCTTGAAGC</th>
      <th>BoneMarrow_5.TGCAATATCTCTACTTAT</th>
      <th>BoneMarrow_5.ACTTATAAAGTTACAATA</th>
      <th>BoneMarrow_5.CACAAGTGCGGATGCGGA</th>
      <th>BoneMarrow_5.CGTATTATTCCATCTACC</th>
      <th>BoneMarrow_5.TGAAGCCCAGACCTCGCA</th>
      <th>BoneMarrow_5.CCTTTCACTTATCTGTGT</th>
      <th>BoneMarrow_5.CGCTTGCACAAGAGGACT</th>
      <th>BoneMarrow_5.TCACTTAGTTTAGCTGTG</th>
      <th>BoneMarrow_5.ACCTGACGTGGCGTATAC</th>
      <th>BoneMarrow_5.TCAAAGCGGCAGGGGCGA</th>
      <th>BoneMarrow_5.ACACCCCTCCATTCTACC</th>
      <th>BoneMarrow_5.GATCTTGGGCGAATGGCG</th>
      <th>BoneMarrow_5.ACTTATGACACTGGCTGC</th>
      <th>BoneMarrow_5.TTTAGGGGCTGCGGCTGC</th>
      <th>BoneMarrow_5.TGATCATCACTTCCAGAC</th>
      <th>BoneMarrow_5.TGCGGAATTCCAGATCTT</th>
      <th>BoneMarrow_5.CAACAAATCTCTCTGAAA</th>
      <th>BoneMarrow_5.ATCAACTTAACTCCAGAC</th>
      <th>BoneMarrow_5.CTCGCATAGAGAGAACGC</th>
      <th>BoneMarrow_5.ACCTGAAACCTAGCGAAT</th>
      <th>BoneMarrow_5.TGATCACCAGACCGAGTA</th>
      <th>BoneMarrow_5.ACGAGCTACTTCACTTAT</th>
      <th>BoneMarrow_5.GATCTTAGGGTCGGACAT</th>
      <th>BoneMarrow_5.GAGGAGGGACATATGCTT</th>
      <th>BoneMarrow_5.TGATCATGATCAAGGACT</th>
      <th>BoneMarrow_5.AAGCGGAACCTAGCTGTG</th>
      <th>BoneMarrow_5.GAGGAGGGACATTCAAAG</th>
      <th>BoneMarrow_5.TGCAATGATCTTGTCCCG</th>
      <th>BoneMarrow_5.ATGGCGAATAAATGCAAT</th>
      <th>BoneMarrow_5.CCTTTCTTTAGGCTGTGT</th>
      <th>BoneMarrow_5.AGGGTCAACCTATGATCA</th>
      <th>BoneMarrow_5.CTGAAAGCAGGATAGAGA</th>
      <th>BoneMarrow_5.CGCTTGCGAGTAGTCGGT</th>
      <th>BoneMarrow_5.CAACAAGTTGCCTATTGT</th>
      <th>BoneMarrow_5.ACCTGACACAAGTTAACT</th>
      <th>BoneMarrow_5.CGTGGCATTTGCGTAATG</th>
      <th>BoneMarrow_5.GCAGGAGTTGCCGGACAT</th>
      <th>BoneMarrow_5.CTCCATTGTCACCGAGTA</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1110008P14Rik</th>
      <td>0.018987</td>
      <td>0.009259</td>
      <td>0.024194</td>
      <td>0.016722</td>
      <td>0.026515</td>
      <td>0.028777</td>
      <td>0.026217</td>
      <td>0.023346</td>
      <td>0.019157</td>
      <td>0.030303</td>
      <td>0.026415</td>
      <td>0.03321</td>
      <td>0.028340</td>
      <td>0.048077</td>
      <td>0.020690</td>
      <td>0.011450</td>
      <td>0.024096</td>
      <td>0.024390</td>
      <td>0.017794</td>
      <td>0.042553</td>
      <td>0.012295</td>
      <td>0.023148</td>
      <td>0.023256</td>
      <td>0.004630</td>
      <td>0.007752</td>
      <td>0.008511</td>
      <td>0.008230</td>
      <td>0.014286</td>
      <td>0.025424</td>
      <td>0.023364</td>
      <td>0.009091</td>
      <td>0.013274</td>
      <td>0.036842</td>
      <td>0.004717</td>
      <td>0.022422</td>
      <td>0.024590</td>
      <td>0.020833</td>
      <td>0.007722</td>
      <td>0.012712</td>
      <td>0.009302</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.25</td>
      <td>0.05</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1700016P03Rik</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1700027J19Rik</th>
      <td>0.012658</td>
      <td>0.009259</td>
      <td>0.008065</td>
      <td>0.020067</td>
      <td>0.015152</td>
      <td>0.007194</td>
      <td>0.011236</td>
      <td>0.000000</td>
      <td>0.003831</td>
      <td>0.000000</td>
      <td>0.003774</td>
      <td>0.00369</td>
      <td>0.012146</td>
      <td>0.004808</td>
      <td>0.010345</td>
      <td>0.007634</td>
      <td>0.004016</td>
      <td>0.009756</td>
      <td>0.014235</td>
      <td>0.015957</td>
      <td>0.016393</td>
      <td>0.013889</td>
      <td>0.009302</td>
      <td>0.023148</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.004115</td>
      <td>0.004762</td>
      <td>0.008475</td>
      <td>0.014019</td>
      <td>0.004545</td>
      <td>0.000000</td>
      <td>0.036842</td>
      <td>0.004717</td>
      <td>0.013453</td>
      <td>0.012295</td>
      <td>0.004167</td>
      <td>0.003861</td>
      <td>0.004237</td>
      <td>0.004651</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1700056N10Rik</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1700119H24Rik</th>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 10000 columns</p>
</div>




```python
mca_counts_norm_lc_output.sum(axis=1)
```




    array([185.41336  ,   5.728099 ,  47.230495 , ...,   5.7145967,
             6.331905 ,   5.6825914], dtype=float32)




```python
mca_counts.head()
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
      <th>Bladder_1.CGGCAGAAAGTTATTCCA</th>
      <th>Bladder_1.CTCGCAAATAAAATCAAC</th>
      <th>Bladder_1.CCATCTAGCGAGTTTAGG</th>
      <th>Bladder_1.GAGGAGCGCTTGATACAG</th>
      <th>Bladder_1.CCAGACACAATAGAATTA</th>
      <th>Bladder_1.CCGACGGGACATATGGCG</th>
      <th>Bladder_1.TAGCATTCAAAGATTCCA</th>
      <th>Bladder_1.CTCCATCCATCTTTTAGG</th>
      <th>Bladder_1.CAAAGTAGGACTAAGTAC</th>
      <th>Bladder_1.TAGCATTGCGGAAACCTA</th>
      <th>Bladder_1.ACACCCACGTTGCACAAG</th>
      <th>Bladder_1.CGCACCATCTCTTAGTCG</th>
      <th>Bladder_1.ACAATAGTCCCGAACGCC</th>
      <th>Bladder_1.CGAGTATAGCATTCAAAG</th>
      <th>Bladder_1.TCGTAATCGGGTGCAGGA</th>
      <th>Bladder_1.ACAATACCGCTATATGTA</th>
      <th>Bladder_1.CGCACCTGATCATTCATA</th>
      <th>Bladder_1.TTCATAGCGAATGTAATG</th>
      <th>Bladder_1.CGCACCTACTTCTCTACC</th>
      <th>Bladder_1.CAACAAGAGATCAACCTA</th>
      <th>Bladder_1.AGGACTCGCACCGAGATC</th>
      <th>Bladder_1.TAGCATCGTATTGCGAAT</th>
      <th>Bladder_1.CCGACGATGCTTCGCTTG</th>
      <th>Bladder_1.TATGTAAAAACGTCTACC</th>
      <th>Bladder_1.CGTATTGATCTTCCTAGA</th>
      <th>Bladder_1.AACCTACAACAAAAAGTT</th>
      <th>Bladder_1.TCACTTGTCGGTTTCATA</th>
      <th>Bladder_1.TTCATATCAAAGTATTGT</th>
      <th>Bladder_1.ACGAGCGAGATCTATTGT</th>
      <th>Bladder_1.AACCTACCATCTGAGGAG</th>
      <th>Bladder_1.AACGCCAAAGTTCCTTTC</th>
      <th>Bladder_1.ACTTATAGTCGTTTTAGG</th>
      <th>Bladder_1.AAGTACATGGCGGCAGGA</th>
      <th>Bladder_1.CATGATAGTCGTCTTCTG</th>
      <th>Bladder_1.AGCGAGAGCGAGAGATGG</th>
      <th>Bladder_1.ACGTTGGTTGCCCATGAT</th>
      <th>Bladder_1.TGATCACACAAGTCAAAG</th>
      <th>Bladder_1.CCTAGATTCCGCGTATAC</th>
      <th>Bladder_1.TTAACTTCGTAAGCAGGA</th>
      <th>Bladder_1.CCAGACATTCCACTGTGT</th>
      <th>...</th>
      <th>BoneMarrow_5.CCTAGACATGATGCTGTG</th>
      <th>BoneMarrow_5.CGCACCATCTCTTGAAGC</th>
      <th>BoneMarrow_5.TGCAATATCTCTACTTAT</th>
      <th>BoneMarrow_5.ACTTATAAAGTTACAATA</th>
      <th>BoneMarrow_5.CACAAGTGCGGATGCGGA</th>
      <th>BoneMarrow_5.CGTATTATTCCATCTACC</th>
      <th>BoneMarrow_5.TGAAGCCCAGACCTCGCA</th>
      <th>BoneMarrow_5.CCTTTCACTTATCTGTGT</th>
      <th>BoneMarrow_5.CGCTTGCACAAGAGGACT</th>
      <th>BoneMarrow_5.TCACTTAGTTTAGCTGTG</th>
      <th>BoneMarrow_5.ACCTGACGTGGCGTATAC</th>
      <th>BoneMarrow_5.TCAAAGCGGCAGGGGCGA</th>
      <th>BoneMarrow_5.ACACCCCTCCATTCTACC</th>
      <th>BoneMarrow_5.GATCTTGGGCGAATGGCG</th>
      <th>BoneMarrow_5.ACTTATGACACTGGCTGC</th>
      <th>BoneMarrow_5.TTTAGGGGCTGCGGCTGC</th>
      <th>BoneMarrow_5.TGATCATCACTTCCAGAC</th>
      <th>BoneMarrow_5.TGCGGAATTCCAGATCTT</th>
      <th>BoneMarrow_5.CAACAAATCTCTCTGAAA</th>
      <th>BoneMarrow_5.ATCAACTTAACTCCAGAC</th>
      <th>BoneMarrow_5.CTCGCATAGAGAGAACGC</th>
      <th>BoneMarrow_5.ACCTGAAACCTAGCGAAT</th>
      <th>BoneMarrow_5.TGATCACCAGACCGAGTA</th>
      <th>BoneMarrow_5.ACGAGCTACTTCACTTAT</th>
      <th>BoneMarrow_5.GATCTTAGGGTCGGACAT</th>
      <th>BoneMarrow_5.GAGGAGGGACATATGCTT</th>
      <th>BoneMarrow_5.TGATCATGATCAAGGACT</th>
      <th>BoneMarrow_5.AAGCGGAACCTAGCTGTG</th>
      <th>BoneMarrow_5.GAGGAGGGACATTCAAAG</th>
      <th>BoneMarrow_5.TGCAATGATCTTGTCCCG</th>
      <th>BoneMarrow_5.ATGGCGAATAAATGCAAT</th>
      <th>BoneMarrow_5.CCTTTCTTTAGGCTGTGT</th>
      <th>BoneMarrow_5.AGGGTCAACCTATGATCA</th>
      <th>BoneMarrow_5.CTGAAAGCAGGATAGAGA</th>
      <th>BoneMarrow_5.CGCTTGCGAGTAGTCGGT</th>
      <th>BoneMarrow_5.CAACAAGTTGCCTATTGT</th>
      <th>BoneMarrow_5.ACCTGACACAAGTTAACT</th>
      <th>BoneMarrow_5.CGTGGCATTTGCGTAATG</th>
      <th>BoneMarrow_5.GCAGGAGTTGCCGGACAT</th>
      <th>BoneMarrow_5.CTCCATTGTCACCGAGTA</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1110008P14Rik</th>
      <td>6</td>
      <td>2</td>
      <td>6</td>
      <td>5</td>
      <td>7</td>
      <td>8</td>
      <td>7</td>
      <td>6</td>
      <td>5</td>
      <td>8</td>
      <td>7</td>
      <td>9</td>
      <td>7</td>
      <td>10</td>
      <td>6</td>
      <td>3</td>
      <td>6</td>
      <td>5</td>
      <td>5</td>
      <td>8</td>
      <td>3</td>
      <td>5</td>
      <td>5</td>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>2</td>
      <td>3</td>
      <td>6</td>
      <td>5</td>
      <td>2</td>
      <td>3</td>
      <td>7</td>
      <td>1</td>
      <td>5</td>
      <td>6</td>
      <td>5</td>
      <td>2</td>
      <td>3</td>
      <td>2</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700016P03Rik</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700027J19Rik</th>
      <td>4</td>
      <td>2</td>
      <td>2</td>
      <td>6</td>
      <td>4</td>
      <td>2</td>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>3</td>
      <td>1</td>
      <td>3</td>
      <td>2</td>
      <td>1</td>
      <td>2</td>
      <td>4</td>
      <td>3</td>
      <td>4</td>
      <td>3</td>
      <td>2</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
      <td>7</td>
      <td>1</td>
      <td>3</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700056N10Rik</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1700119H24Rik</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 10000 columns</p>
</div>








```python

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
      <th>GeneName</th>
      <th>Bladder_1.CGGCAGAAAGTTATTCCA</th>
      <th>Bladder_1.CTCGCAAATAAAATCAAC</th>
      <th>Bladder_1.CCATCTAGCGAGTTTAGG</th>
      <th>Bladder_1.GAGGAGCGCTTGATACAG</th>
      <th>Bladder_1.CCAGACACAATAGAATTA</th>
      <th>Bladder_1.CCGACGGGACATATGGCG</th>
      <th>Bladder_1.TAGCATTCAAAGATTCCA</th>
      <th>Bladder_1.CTCCATCCATCTTTTAGG</th>
      <th>Bladder_1.CAAAGTAGGACTAAGTAC</th>
      <th>Bladder_1.TAGCATTGCGGAAACCTA</th>
      <th>Bladder_1.ACACCCACGTTGCACAAG</th>
      <th>Bladder_1.CGCACCATCTCTTAGTCG</th>
      <th>Bladder_1.ACAATAGTCCCGAACGCC</th>
      <th>Bladder_1.CGAGTATAGCATTCAAAG</th>
      <th>Bladder_1.TCGTAATCGGGTGCAGGA</th>
      <th>Bladder_1.ACAATACCGCTATATGTA</th>
      <th>Bladder_1.CGCACCTGATCATTCATA</th>
      <th>Bladder_1.TTCATAGCGAATGTAATG</th>
      <th>Bladder_1.CGCACCTACTTCTCTACC</th>
      <th>Bladder_1.CAACAAGAGATCAACCTA</th>
      <th>Bladder_1.AGGACTCGCACCGAGATC</th>
      <th>Bladder_1.TAGCATCGTATTGCGAAT</th>
      <th>Bladder_1.CCGACGATGCTTCGCTTG</th>
      <th>Bladder_1.TATGTAAAAACGTCTACC</th>
      <th>Bladder_1.CGTATTGATCTTCCTAGA</th>
      <th>Bladder_1.AACCTACAACAAAAAGTT</th>
      <th>Bladder_1.TCACTTGTCGGTTTCATA</th>
      <th>Bladder_1.TTCATATCAAAGTATTGT</th>
      <th>Bladder_1.ACGAGCGAGATCTATTGT</th>
      <th>Bladder_1.AACCTACCATCTGAGGAG</th>
      <th>Bladder_1.AACGCCAAAGTTCCTTTC</th>
      <th>Bladder_1.ACTTATAGTCGTTTTAGG</th>
      <th>Bladder_1.AAGTACATGGCGGCAGGA</th>
      <th>Bladder_1.CATGATAGTCGTCTTCTG</th>
      <th>Bladder_1.AGCGAGAGCGAGAGATGG</th>
      <th>Bladder_1.ACGTTGGTTGCCCATGAT</th>
      <th>Bladder_1.TGATCACACAAGTCAAAG</th>
      <th>Bladder_1.CCTAGATTCCGCGTATAC</th>
      <th>Bladder_1.TTAACTTCGTAAGCAGGA</th>
      <th>...</th>
      <th>BoneMarrow_5.CCTAGACATGATGCTGTG</th>
      <th>BoneMarrow_5.CGCACCATCTCTTGAAGC</th>
      <th>BoneMarrow_5.TGCAATATCTCTACTTAT</th>
      <th>BoneMarrow_5.ACTTATAAAGTTACAATA</th>
      <th>BoneMarrow_5.CACAAGTGCGGATGCGGA</th>
      <th>BoneMarrow_5.CGTATTATTCCATCTACC</th>
      <th>BoneMarrow_5.TGAAGCCCAGACCTCGCA</th>
      <th>BoneMarrow_5.CCTTTCACTTATCTGTGT</th>
      <th>BoneMarrow_5.CGCTTGCACAAGAGGACT</th>
      <th>BoneMarrow_5.TCACTTAGTTTAGCTGTG</th>
      <th>BoneMarrow_5.ACCTGACGTGGCGTATAC</th>
      <th>BoneMarrow_5.TCAAAGCGGCAGGGGCGA</th>
      <th>BoneMarrow_5.ACACCCCTCCATTCTACC</th>
      <th>BoneMarrow_5.GATCTTGGGCGAATGGCG</th>
      <th>BoneMarrow_5.ACTTATGACACTGGCTGC</th>
      <th>BoneMarrow_5.TTTAGGGGCTGCGGCTGC</th>
      <th>BoneMarrow_5.TGATCATCACTTCCAGAC</th>
      <th>BoneMarrow_5.TGCGGAATTCCAGATCTT</th>
      <th>BoneMarrow_5.CAACAAATCTCTCTGAAA</th>
      <th>BoneMarrow_5.ATCAACTTAACTCCAGAC</th>
      <th>BoneMarrow_5.CTCGCATAGAGAGAACGC</th>
      <th>BoneMarrow_5.ACCTGAAACCTAGCGAAT</th>
      <th>BoneMarrow_5.TGATCACCAGACCGAGTA</th>
      <th>BoneMarrow_5.ACGAGCTACTTCACTTAT</th>
      <th>BoneMarrow_5.GATCTTAGGGTCGGACAT</th>
      <th>BoneMarrow_5.GAGGAGGGACATATGCTT</th>
      <th>BoneMarrow_5.TGATCATGATCAAGGACT</th>
      <th>BoneMarrow_5.AAGCGGAACCTAGCTGTG</th>
      <th>BoneMarrow_5.GAGGAGGGACATTCAAAG</th>
      <th>BoneMarrow_5.TGCAATGATCTTGTCCCG</th>
      <th>BoneMarrow_5.ATGGCGAATAAATGCAAT</th>
      <th>BoneMarrow_5.CCTTTCTTTAGGCTGTGT</th>
      <th>BoneMarrow_5.AGGGTCAACCTATGATCA</th>
      <th>BoneMarrow_5.CTGAAAGCAGGATAGAGA</th>
      <th>BoneMarrow_5.CGCTTGCGAGTAGTCGGT</th>
      <th>BoneMarrow_5.CAACAAGTTGCCTATTGT</th>
      <th>BoneMarrow_5.ACCTGACACAAGTTAACT</th>
      <th>BoneMarrow_5.CGTGGCATTTGCGTAATG</th>
      <th>BoneMarrow_5.GCAGGAGTTGCCGGACAT</th>
      <th>BoneMarrow_5.CTCCATTGTCACCGAGTA</th>
    </tr>
    <tr>
      <th>GeneName</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1110008P14Rik</td>
      <td>0.018987</td>
      <td>0.009259</td>
      <td>0.024194</td>
      <td>0.016722</td>
      <td>0.026515</td>
      <td>0.028777</td>
      <td>0.026217</td>
      <td>0.023346</td>
      <td>0.019157</td>
      <td>0.030303</td>
      <td>0.026415</td>
      <td>0.03321</td>
      <td>0.028340</td>
      <td>0.048077</td>
      <td>0.020690</td>
      <td>0.011450</td>
      <td>0.024096</td>
      <td>0.024390</td>
      <td>0.017794</td>
      <td>0.042553</td>
      <td>0.012295</td>
      <td>0.023148</td>
      <td>0.023256</td>
      <td>0.004630</td>
      <td>0.007752</td>
      <td>0.008511</td>
      <td>0.008230</td>
      <td>0.014286</td>
      <td>0.025424</td>
      <td>0.023364</td>
      <td>0.009091</td>
      <td>0.013274</td>
      <td>0.036842</td>
      <td>0.004717</td>
      <td>0.022422</td>
      <td>0.024590</td>
      <td>0.020833</td>
      <td>0.007722</td>
      <td>0.012712</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.25</td>
      <td>0.05</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1700016P03Rik</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1700027J19Rik</td>
      <td>0.012658</td>
      <td>0.009259</td>
      <td>0.008065</td>
      <td>0.020067</td>
      <td>0.015152</td>
      <td>0.007194</td>
      <td>0.011236</td>
      <td>0.000000</td>
      <td>0.003831</td>
      <td>0.000000</td>
      <td>0.003774</td>
      <td>0.00369</td>
      <td>0.012146</td>
      <td>0.004808</td>
      <td>0.010345</td>
      <td>0.007634</td>
      <td>0.004016</td>
      <td>0.009756</td>
      <td>0.014235</td>
      <td>0.015957</td>
      <td>0.016393</td>
      <td>0.013889</td>
      <td>0.009302</td>
      <td>0.023148</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.004115</td>
      <td>0.004762</td>
      <td>0.008475</td>
      <td>0.014019</td>
      <td>0.004545</td>
      <td>0.000000</td>
      <td>0.036842</td>
      <td>0.004717</td>
      <td>0.013453</td>
      <td>0.012295</td>
      <td>0.004167</td>
      <td>0.003861</td>
      <td>0.004237</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1700056N10Rik</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1700119H24Rik</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1995</th>
      <td>Slc22a3</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1996</th>
      <td>Treml1</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1997</th>
      <td>Adamts13</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1998</th>
      <td>C8a</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1999</th>
      <td>Ighv5-17</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.00000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>...</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
<p>2000 rows × 10001 columns</p>
</div>




```python
from google.colab import files
mca_counts_norm_lc.to_csv('mca_counts_norm_lc.csv')
files.download('mca_counts_norm_lc.csv')
```


```python
mca_counts_norm_lc_output_pd = pd.DataFrame(data=mca_counts_norm_lc_output[0:,0:])  # # # values 1st column as index 1st row as the column names
mca_counts_norm_lc_output_pd.columns = mca_counts_norm_lc.columns[1:]

from google.colab import files
mca_counts_norm_lc_output_pd.to_csv('mca_counts_norm_lc_output_pd_colnames.csv')
files.download('mca_counts_norm_lc_output_pd_colnames.csv')
```


```r
mca.10K.original <- read.csv("C:\\Sivome\\SingleCellGenomics\\MCA\\MCA\\mca_counts_norm_lc.csv")
rownames(mca.10K.original) <- mca.10K.original$X
mca.10K.original <- mca.10K.original[, -1]


mca.10K.autoencoder <- read.csv("C:\\Sivome\\SingleCellGenomics\\MCA\\MCA\\mca_counts_norm_lc_output_pd_colnames.csv")
mca.10K.autoencoder <- mca.10K.autoencoder[, -1]
rownames(mca.10K.autoencoder) <- rownames(mca.10K.original)

mca.10K.original <- data.matrix(mca.10K.original)
mca.10K.autoencoder <- data.matrix(mca.10K.autoencoder)
```


```r
# duplicate mca.10K seurat objects and replace the scale.data with the new original and autoencoder analyses
mca.10K.original.seurat <- mca.10K
mca.10K.autoencoder.seurat <- mca.10K
mca.10K.original.seurat@assays$RNA@scale.data <- mca.10K.original
mca.10K.autoencoder.seurat@assays$RNA@scale.data <- mca.10K.autoencoder
```


```r
mca.10K.original.seurat <- RunPCA(mca.10K.original.seurat, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
```

```
## PC_ 1
## Positive:  Ighv2-2, Ighv1-22, Ighv10-3, Ighv3-6, Ighv1-81
## Negative:  Retnlg, S100a6, mt-Cytb, Rps27, Rplp1
## PC_ 2
## Positive:  Retnlg, Elane, Mpo, Fcnb, Igkc
## Negative:  Mgp, Dcn, Gstm1, Gsn, S100a6
## PC_ 3
## Positive:  Retnlg, Mgp, Dcn, Gsn, Col1a2
## Negative:  Gstm1, Elane, Rplp0, Rplp1, mt-Cytb
## PC_ 4
## Positive:  Retnlg, S100a6, Gstm1, Gsta4, Ly6d
## Negative:  Elane, Mpo, Mgp, Dcn, Gsn
## PC_ 5
## Positive:  Crip1, Rplp0, Rplp1, Rps27, Rps29
## Negative:  Elane, Gstm1, S100a6, Mpo, Prtn3
```

```r
mca.10K.original.seurat <- FindNeighbors(mca.10K.original.seurat, reduction = "pca", dims = 1:75, nn.eps = 0.5)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
mca.10K.original.seurat <- FindClusters(mca.10K.original.seurat, resolution = 3, n.start = 10)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
##
## Number of nodes: 10000
## Number of edges: 486341
##
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.7731
## Number of communities: 36
## Elapsed time: 1 seconds
```

```r
#t-SNE, a populare classification method
mca.10K.original.seurat <- RunTSNE(mca.10K.original.seurat, dims = 1:75)

#UMAP is relatively new and with some datasets, it is shown to perform better than t-SNE
mca.10K.original.seurat <- RunUMAP(mca.10K.original.seurat, dims = 1:75, min.dist = 0.75)
```

```
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
```

```
## 20:07:02 UMAP embedding parameters a = 0.2734 b = 1.622
```

```
## 20:07:02 Read 10000 rows and found 75 numeric columns
```

```
## 20:07:02 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 20:07:02 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 20:07:05 Writing NN index file to temp file C:\Users\Viswa\AppData\Local\Temp\Rtmpm6apjQ\filee7043871a7
## 20:07:05 Searching Annoy index using 1 thread, search_k = 3000
## 20:07:08 Annoy recall = 100%
## 20:07:08 Commencing smooth kNN distance calibration using 1 thread
## 20:07:09 Initializing from normalized Laplacian + noise
## 20:07:11 Commencing optimization for 500 epochs, with 482572 positive edges
## 20:07:36 Optimization finished
```


```r
> DimPlot(mca.10K.original.seurat, group.by = 'orig.ident') + DarkTheme()
> DimPlot(mca.10K.autoencoder.seurat, group.by = 'orig.ident') + DarkTheme()
```

```
## Error: <text>:1:1: unexpected '>'
## 1: >
##     ^
```




```r
mca.10K.autoencoder.seurat <- RunPCA(mca.10K.autoencoder.seurat, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
```

```
## PC_ 1
## Positive:  Gm12583, Llph-ps1, Gm5356, Arg1, Gm4925
## Negative:  Mgp, Gsn, Sparc, Col1a2, Col3a1
## PC_ 2
## Positive:  Mgp, Gsn, Col1a2, Sparc, Col3a1
## Negative:  Plac8, Rplp0, mt-Cytb, Rplp1, Mpo
## PC_ 3
## Positive:  Gsta4, Sprr1a, Krt15, Wfdc2, Ly6d
## Negative:  Dcn, Col1a2, Col3a1, Mgp, Elane
## PC_ 4
## Positive:  Gsta4, Sprr1a, Krt15, Mgp, Wfdc2
## Negative:  Acta2, Myl9, Igfbp7, Tpm2, Hspb1
## PC_ 5
## Positive:  Elane, Plac8, Rplp0, Mpo, mt-Cytb
## Negative:  Hbb-bt, Hba-a1, Car2, Car1, Blvrb
```

```r
mca.10K.autoencoder.seurat <- FindNeighbors(mca.10K.autoencoder.seurat, reduction = "pca", dims = 1:75, nn.eps = 0.5)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
mca.10K.autoencoder.seurat <- FindClusters(mca.10K.autoencoder.seurat, resolution = 3, n.start = 10)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
##
## Number of nodes: 10000
## Number of edges: 300384
##
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8062
## Number of communities: 41
## Elapsed time: 0 seconds
```

```r
#t-SNE, a populare classification method
mca.10K.autoencoder.seurat <- RunTSNE(mca.10K.autoencoder.seurat, dims = 1:75)

#UMAP is relatively new and with some datasets, it is shown to perform better than t-SNE
mca.10K.autoencoder.seurat <- RunUMAP(mca.10K.autoencoder.seurat, dims = 1:75, min.dist = 0.75)
```

```
## 20:08:12 UMAP embedding parameters a = 0.2734 b = 1.622
```

```
## 20:08:12 Read 10000 rows and found 75 numeric columns
```

```
## 20:08:12 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 20:08:12 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 20:08:14 Writing NN index file to temp file C:\Users\Viswa\AppData\Local\Temp\Rtmpm6apjQ\filee70577043cf
## 20:08:14 Searching Annoy index using 1 thread, search_k = 3000
## 20:08:16 Annoy recall = 100%
## 20:08:17 Commencing smooth kNN distance calibration using 1 thread
## 20:08:18 Initializing from normalized Laplacian + noise
## 20:08:18 Commencing optimization for 500 epochs, with 395680 positive edges
## 20:08:41 Optimization finished
```
