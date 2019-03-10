---
title: "Reporting Proteomics Results with Jupyter Notebook"
date: '2019-03-08'
layout: post
categories: Proteomics
---

As a continuation to my previous blog on Proteomics with OMSSA, here I extend the analyses to look at a open-access yeast data set. I used only one raw file in the entire project and filtered it further to use few thousands of spectra. This is for blog, so I have the benefit to minimize the data for faster processing.

Few open-source efforts to get the mass-spectrometry based proteomics efforts are [Proteome Exchange Consortium], [PRIDE]. I used trial version of a commercial software to convert it for a format this is very easy to analyze. The format is [MGF].

Before we dive into the analyses, we need few things.
1. Fasta file for the species that we are analyzing. Since this is Yeast data, we need Yeast fasta file. There are multiple places to get this fasta and in the previous blogs, I mentioned about the NCBI resources. This time, let's focus on [UniProt](https://www.uniprot.org/). This is a resource, as the name says, provides information on proteins that is comprehensive, well annotated. There are different versions of the database depending on the kinds of analyses you're doing. Here, I used [Baker's yeast](https://www.uniprot.org/taxonomy/559292). There are few options of which variant of the database you can pick in UniProt. One such variant is SwissProt that has reviewed protein sequences. Here is how the fasta file looks like:

```console
C:\Users\Viswa\blastDb>head -n 20 S288c.fasta
>sp|P32768|FLO1_YEAST Flocculation protein FLO1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=FLO1 PE=1 SV=4
MTMPHRYMFLAVFTLLALTSVASGATEACLPAGQRKSGMNINFYQYSLKDSSTYSNAAYM
AYGYASKTKLGSVGGQTDISIDYNIPCVSSSGTFPCPQEDSYGNWGCKGMGACSNSQGIA
YWSTDLFGFYTTPTNVTLEMTGYFLPPQTGSYTFKFATVDDSAILSVGGATAFNCCAQQQ
PPITSTNFTIDGIKPWGGSLPPNIEGTVYMYAGYYYPMKVVYSNAVSWGTLPISVTLPDG
TTVSDDFEGYVYSFDDDLSQSNCTVPDPSNYAVSTTTTTTEPWTGTFTSTSTEMTTVTGT
NGVPTDETVIVIRTPTTASTIITTTEPWNSTFTSTSTELTTVTGTNGVRTDETIIVIRTP
TTATTAITTTEPWNSTFTSTSTELTTVTGTNGLPTDETIIVIRTPTTATTAMTTTQPWND
TFTSTSTELTTVTGTNGLPTDETIIVIRTPTTATTAMTTTQPWNDTFTSTSTELTTVTGT
NGLPTDETIIVIRTPTTATTAMTTTQPWNDTFTSTSTEITTVTGTNGLPTDETIIVIRTP
TTATTAMTTPQPWNDTFTSTSTEMTTVTGTNGLPTDETIIVIRTPTTATTAITTTEPWNS
TFTSTSTEMTTVTGTNGLPTDETIIVIRTPTTATTAITTTQPWNDTFTSTSTEMTTVTGT
NGLPTDETIIVIRTPTTATTAMTTTQPWNDTFTSTSTEITTVTGTTGLPTDETIIVIRTP
TTATTAMTTTQPWNDTFTSTSTEMTTVTGTNGVPTDETVIVIRTPTSEGLISTTTEPWTG
TFTSTSTEMTTVTGTNGQPTDETVIVIRTPTSEGLVTTTTEPWTGTFTSTSTEMTTITGT
NGVPTDETVIVIRTPTSEGLISTTTEPWTGTFTSTSTEMTTITGTNGQPTDETVIVIRTP
TSEGLISTTTEPWTGTFTSTSTEMTHVTGTNGVPTDETVIVIRTPTSEGLISTTTEPWTG
TFTSTSTEVTTITGTNGQPTDETVIVIRTPTSEGLISTTTEPWTGTFTSTSTEMTTVTGT
NGQPTDETVIVIRTPTSEGLVTTTTEPWTGTFTSTSTEMSTVTGTNGLPTDETVIVVKTP
TTAISSSLSSSSSGQITSSITSSRPIITPFYPSNGTSVISSSVISSSVTSSLFTSSPVIS

C:\Users\Viswa\blastDb>
```

You can download the fasta file [here](https://www.uniprot.org/uniprot/?query=S288c+AND+reviewed%3Ayes&sort=score).
if you read the [earlier blog](https://sivome.github.io/proteomics/2019/03/02/Proteomics-with-OMSSA.html), it is easy to figure out how to use the fasta file for proteomics search with OMSSA. There is something called target-decoy database that you need to use for proteomics search. As the name says, we need to append a decoy version to the target fasta file so that search algorithm can use the metrics from decoy hits for down-stream analyses. More information on the target-decoy database search [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922680/). [There is a open-source perl script](https://edwardslab.bmcb.georgetown.edu/mascot/help/decoy_help.html) available that does this job of creating decoy version of the fasta file and appends the decoy to the target fasta file. The decoy version, in the current blog, will have ###REV### appended to the protein name.

Now that we have the database available to see which parts of the protein can be mapped with the mass-spec data, let's take a look at the actual data from the mass-spec. As I mentioned earlier, the current dataset used is MGF version i.e., mascot generic version. This is how a sample MGF file looks like:

```console
C:\Users\Viswa\folder_for_analyses>head -n 20 Yeast.mgf
BEGIN IONS
TITLE=index=0
PEPMASS=519.7517
CHARGE=2+
SCANS=16
RTINSECONDS=1.878
149.12 30.216585
149.985 29.74796
175.185 20.990004
177.198 17.410294
201.096 16.743763
230.869 28.351389
243.217 27.75121
255.166 10.026513
272.247 34.915855
273.127 19.958282
281.883 17.084764
290.254 20.235743
294.149 14.525609
300.068 33.84702
```
As you can see, the peptide mass [PEPMASS], charge state [CHARGE], scan number, retention time and the peak lists are present in this format. OMSSA uses the PEPMASS*CHARGE as possible observed peptide mass, and then looks for the peptide that matches this mass. In this case, the enzyme used to digest the protein is trypsin. Since trypsin cleaves at K/R, the program digests the protein sequence (both for targets and decoys), and if the in-silico peptide from the digest matches the PEPMASS*CHARGE from MGF file, it scores the peptide. In OMSSA, you can think of this score as E-value (similar to BLAST scoring). It does this for all the scans in the mgf file and scores all the possible peptide matches. This can be understood easily from the analyses below.

Let's skip the actual proteomics search and focus on the output of the OMSSA results. Let's say we have this in csv format and named it S288c_run.csv.

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
```


```python
# Display multiple columns, common in high-dimensional space (add rows as well)
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 20)
```

Here, we read the OMSSA output into a pandas dataframe.

```python
omssa_output = pd.read_csv("S288c_run.csv")
```

Let's see how this looks like:
```python
omssa_output.head() # Talk about different accessions?
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
      <th>Spectrum number</th>
      <th>Filename/id</th>
      <th>Peptide</th>
      <th>E-value</th>
      <th>Mass</th>
      <th>gi</th>
      <th>Accession</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Defline</th>
      <th>Mods</th>
      <th>Charge</th>
      <th>Theo Mass</th>
      <th>P-value</th>
      <th>NIST score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>120</td>
      <td>index=120</td>
      <td>GSIDEQHPR</td>
      <td>0.099394</td>
      <td>1037.487</td>
      <td>0</td>
      <td>BL_ORD_ID:50</td>
      <td>250</td>
      <td>258</td>
      <td>sp|P06169|PDC1_YEAST Pyruvate decarboxylase is...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1037.490</td>
      <td>3.803824e-05</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>120</td>
      <td>index=120</td>
      <td>GSIDEQHPR</td>
      <td>0.099394</td>
      <td>1037.487</td>
      <td>0</td>
      <td>BL_ORD_ID:131</td>
      <td>250</td>
      <td>258</td>
      <td>sp|P26263|PDC6_YEAST Pyruvate decarboxylase is...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1037.490</td>
      <td>3.803824e-05</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>339</td>
      <td>index=339</td>
      <td>TSGRPIKGDSSAGGK</td>
      <td>0.001174</td>
      <td>1416.728</td>
      <td>0</td>
      <td>BL_ORD_ID:2647</td>
      <td>176</td>
      <td>190</td>
      <td>sp|P47075|VTC4_YEAST Vacuolar transporter chap...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1416.731</td>
      <td>5.359808e-07</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>421</td>
      <td>index=421</td>
      <td>TSGRPIKGDSSAGGK</td>
      <td>0.035522</td>
      <td>1416.729</td>
      <td>0</td>
      <td>BL_ORD_ID:2647</td>
      <td>176</td>
      <td>190</td>
      <td>sp|P47075|VTC4_YEAST Vacuolar transporter chap...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1416.731</td>
      <td>1.624232e-05</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>506</td>
      <td>index=506</td>
      <td>NEETSGEGGEDKNEPSSK</td>
      <td>0.028961</td>
      <td>1892.786</td>
      <td>0</td>
      <td>BL_ORD_ID:4527</td>
      <td>76</td>
      <td>93</td>
      <td>sp|Q02776|TIM50_YEAST Mitochondrial import inn...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1892.789</td>
      <td>1.991815e-05</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
omssa_output.shape
```
It looks like there are a bit more than 1700 peptide spectral matches [PSMs].



    (1703, 16)




```python
omssa_output.columns
```

These are the columns of the OMSSA output.


    Index(['Spectrum number', ' Filename/id', ' Peptide', ' E-value', ' Mass',
           ' gi', ' Accession', ' Start', ' Stop', ' Defline', ' Mods', ' Charge',
           ' Theo Mass', ' P-value', ' NIST score'],
          dtype='object')


Let's take a look at the E-value distribution.

```python
eps=1e-32 # To take care of zeoes
fig, ax = plt.subplots()
plt.hist(np.log10(eps + omssa_output[' E-value']), bins = 100)
loc = ticker.MultipleLocator(base=2) # thanks again to google!
ax.xaxis.set_major_locator(loc)
ax.set_xlabel('Score in E-value')
ax.set_ylabel('Counts')
ax.set_title('Score distribution')
fig.tight_layout
```




    <bound method Figure.tight_layout of <Figure size 432x288 with 1 Axes>>




![png](output_6_1.png)

We searched for oxidation of methionine as a variable modification. Let's see how many peptides are oxidized with methionine.

```python
# Look at oxidation
omssa_output[' IsMod'] = omssa_output[' Mods'].str.match('oxidation', na=False)
```


```python
mods_only = omssa_output[omssa_output[' IsMod'].values]
```

Here is the table of the hits with modifications.
```python
# Since there are few modifications, display all # Note the REV on one of the hits, this is coming from the reverse database.
mods_only
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
      <th>Spectrum number</th>
      <th>Filename/id</th>
      <th>Peptide</th>
      <th>E-value</th>
      <th>Mass</th>
      <th>gi</th>
      <th>Accession</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Defline</th>
      <th>Mods</th>
      <th>Charge</th>
      <th>Theo Mass</th>
      <th>P-value</th>
      <th>NIST score</th>
      <th>IsReverse</th>
      <th>IsMod</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>55</th>
      <td>729</td>
      <td>index=729</td>
      <td>SKQEASQmAAmAEK</td>
      <td>3.606349e-03</td>
      <td>1540.685</td>
      <td>0</td>
      <td>BL_ORD_ID:6013</td>
      <td>655</td>
      <td>668</td>
      <td>sp|P32589|HSP7F_YEAST Heat shock protein homol...</td>
      <td>oxidation of M:8 ,oxidation of M:11</td>
      <td>2</td>
      <td>1540.687</td>
      <td>2.128896e-06</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>91</th>
      <td>1001</td>
      <td>index=1001</td>
      <td>YATmTGHHVER</td>
      <td>2.314726e-02</td>
      <td>1316.591</td>
      <td>0</td>
      <td>BL_ORD_ID:4636</td>
      <td>70</td>
      <td>80</td>
      <td>sp|P09436|SYIC_YEAST Isoleucine--tRNA ligase, ...</td>
      <td>oxidation of M:4</td>
      <td>3</td>
      <td>1316.593</td>
      <td>1.120390e-05</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>92</th>
      <td>1005</td>
      <td>index=1005</td>
      <td>YATmTGHHVER</td>
      <td>6.635733e-03</td>
      <td>1316.593</td>
      <td>0</td>
      <td>BL_ORD_ID:4636</td>
      <td>70</td>
      <td>80</td>
      <td>sp|P09436|SYIC_YEAST Isoleucine--tRNA ligase, ...</td>
      <td>oxidation of M:4</td>
      <td>2</td>
      <td>1316.593</td>
      <td>3.218105e-06</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>243</th>
      <td>1754</td>
      <td>index=1754</td>
      <td>TPAEmSRPATTTR</td>
      <td>4.902549e-01</td>
      <td>1433.174</td>
      <td>0</td>
      <td>BL_ORD_ID:3809</td>
      <td>1361</td>
      <td>1373</td>
      <td>sp|P19097|FAS2_YEAST Fatty acid synthase subun...</td>
      <td>oxidation of M:5</td>
      <td>3</td>
      <td>1433.695</td>
      <td>2.362674e-04</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>255</th>
      <td>1823</td>
      <td>index=1823</td>
      <td>TPAEmSRPATTTR</td>
      <td>5.176514e-01</td>
      <td>1433.690</td>
      <td>0</td>
      <td>BL_ORD_ID:3809</td>
      <td>1361</td>
      <td>1373</td>
      <td>sp|P19097|FAS2_YEAST Fatty acid synthase subun...</td>
      <td>oxidation of M:5</td>
      <td>3</td>
      <td>1433.695</td>
      <td>2.499524e-04</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>352</th>
      <td>2220</td>
      <td>index=2220</td>
      <td>AITmTNTNEVAGNTTCNKSR</td>
      <td>8.615245e-01</td>
      <td>2198.522</td>
      <td>0</td>
      <td>BL_ORD_ID:396</td>
      <td>640</td>
      <td>659</td>
      <td>sp|P38148|PPS1_YEAST Dual specificity protein ...</td>
      <td>oxidation of M:4</td>
      <td>3</td>
      <td>2198.006</td>
      <td>6.443713e-04</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>543</th>
      <td>2907</td>
      <td>index=2907</td>
      <td>QINENDAEAmNK</td>
      <td>7.328834e-02</td>
      <td>1391.597</td>
      <td>0</td>
      <td>BL_ORD_ID:3789</td>
      <td>775</td>
      <td>786</td>
      <td>sp|P53978|EF3B_YEAST Elongation factor 3B OS=S...</td>
      <td>oxidation of M:10</td>
      <td>2</td>
      <td>1391.600</td>
      <td>3.612042e-05</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>544</th>
      <td>2907</td>
      <td>index=2907</td>
      <td>QINENDAEAmNK</td>
      <td>7.328834e-02</td>
      <td>1391.597</td>
      <td>0</td>
      <td>BL_ORD_ID:5834</td>
      <td>775</td>
      <td>786</td>
      <td>sp|P16521|EF3A_YEAST Elongation factor 3A OS=S...</td>
      <td>oxidation of M:10</td>
      <td>2</td>
      <td>1391.600</td>
      <td>3.612042e-05</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>545</th>
      <td>2909</td>
      <td>index=2909</td>
      <td>SKQEASQmAAMAEK</td>
      <td>3.483283e-07</td>
      <td>1524.693</td>
      <td>0</td>
      <td>BL_ORD_ID:6013</td>
      <td>655</td>
      <td>668</td>
      <td>sp|P32589|HSP7F_YEAST Heat shock protein homol...</td>
      <td>oxidation of M:8</td>
      <td>2</td>
      <td>1524.692</td>
      <td>1.949235e-10</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>546</th>
      <td>2909</td>
      <td>index=2909</td>
      <td>SKQEASQMAAmAEK</td>
      <td>6.834636e-03</td>
      <td>1524.693</td>
      <td>0</td>
      <td>BL_ORD_ID:6013</td>
      <td>655</td>
      <td>668</td>
      <td>sp|P32589|HSP7F_YEAST Heat shock protein homol...</td>
      <td>oxidation of M:11</td>
      <td>2</td>
      <td>1524.692</td>
      <td>3.824642e-06</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
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
    </tr>
    <tr>
      <th>901</th>
      <td>5994</td>
      <td>index=5994</td>
      <td>TRLGLLIGmIIGVQK</td>
      <td>6.219942e-01</td>
      <td>1628.717</td>
      <td>0</td>
      <td>BL_ORD_ID:8146</td>
      <td>65</td>
      <td>79</td>
      <td>###REV###sp|P47118|YAE1_YEAST Reverse sequence...</td>
      <td>oxidation of M:9</td>
      <td>3</td>
      <td>1626.984</td>
      <td>3.408187e-04</td>
      <td>0</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>956</th>
      <td>6357</td>
      <td>index=6357</td>
      <td>mGHSGAIVEGSGTDAESKK</td>
      <td>1.819108e-03</td>
      <td>1875.860</td>
      <td>0</td>
      <td>BL_ORD_ID:6570</td>
      <td>282</td>
      <td>300</td>
      <td>sp|P53598|SUCA_YEAST Succinate--CoA ligase [AD...</td>
      <td>oxidation of M:1</td>
      <td>3</td>
      <td>1875.862</td>
      <td>1.202319e-06</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1006</th>
      <td>6733</td>
      <td>index=6733</td>
      <td>HQGImVGmGQK</td>
      <td>1.049045e-03</td>
      <td>1216.567</td>
      <td>0</td>
      <td>BL_ORD_ID:4846</td>
      <td>40</td>
      <td>50</td>
      <td>sp|P60010|ACT_YEAST Actin OS=Saccharomyces cer...</td>
      <td>oxidation of M:5 ,oxidation of M:8</td>
      <td>2</td>
      <td>1216.568</td>
      <td>4.465924e-07</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1184</th>
      <td>7854</td>
      <td>index=7854</td>
      <td>HVGDmGNVKTDENGVAK</td>
      <td>1.617043e-04</td>
      <td>1785.827</td>
      <td>0</td>
      <td>BL_ORD_ID:4184</td>
      <td>81</td>
      <td>97</td>
      <td>sp|P00445|SODC_YEAST Superoxide dismutase [Cu-...</td>
      <td>oxidation of M:5</td>
      <td>3</td>
      <td>1785.830</td>
      <td>1.029964e-07</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1185</th>
      <td>7855</td>
      <td>index=7855</td>
      <td>HVGDmGNVKTDENGVAK</td>
      <td>4.194812e-03</td>
      <td>1785.827</td>
      <td>0</td>
      <td>BL_ORD_ID:4184</td>
      <td>81</td>
      <td>97</td>
      <td>sp|P00445|SODC_YEAST Superoxide dismutase [Cu-...</td>
      <td>oxidation of M:5</td>
      <td>3</td>
      <td>1785.830</td>
      <td>2.671854e-06</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1334</th>
      <td>8769</td>
      <td>index=8769</td>
      <td>DLmACAQTGSGK</td>
      <td>4.755882e-01</td>
      <td>1253.535</td>
      <td>0</td>
      <td>BL_ORD_ID:3958</td>
      <td>193</td>
      <td>204</td>
      <td>sp|P24784|DBP1_YEAST ATP-dependent RNA helicas...</td>
      <td>oxidation of M:3</td>
      <td>2</td>
      <td>1253.537</td>
      <td>2.295310e-04</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1335</th>
      <td>8769</td>
      <td>index=8769</td>
      <td>DLmACAQTGSGK</td>
      <td>4.755882e-01</td>
      <td>1253.535</td>
      <td>0</td>
      <td>BL_ORD_ID:4004</td>
      <td>181</td>
      <td>192</td>
      <td>sp|P06634|DED1_YEAST ATP-dependent RNA helicas...</td>
      <td>oxidation of M:3</td>
      <td>2</td>
      <td>1253.537</td>
      <td>2.295310e-04</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1336</th>
      <td>8774</td>
      <td>index=8774</td>
      <td>SKQEASQMAAmAEK</td>
      <td>1.248135e-04</td>
      <td>1524.689</td>
      <td>0</td>
      <td>BL_ORD_ID:6013</td>
      <td>655</td>
      <td>668</td>
      <td>sp|P32589|HSP7F_YEAST Heat shock protein homol...</td>
      <td>oxidation of M:11</td>
      <td>2</td>
      <td>1524.692</td>
      <td>7.023830e-08</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1493</th>
      <td>9350</td>
      <td>index=9350</td>
      <td>TLTAQSmQNSTQSAPNK</td>
      <td>5.296256e-05</td>
      <td>1821.851</td>
      <td>0</td>
      <td>BL_ORD_ID:133</td>
      <td>49</td>
      <td>65</td>
      <td>sp|P33302|PDR5_YEAST Pleiotropic ABC efflux tr...</td>
      <td>oxidation of M:7</td>
      <td>2</td>
      <td>1821.855</td>
      <td>3.213748e-08</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1675</th>
      <td>9904</td>
      <td>index=9904</td>
      <td>TQEFSQmACEK</td>
      <td>1.701215e-01</td>
      <td>1373.555</td>
      <td>0</td>
      <td>BL_ORD_ID:4512</td>
      <td>186</td>
      <td>196</td>
      <td>sp|Q04947|RTN1_YEAST Reticulon-like protein 1 ...</td>
      <td>oxidation of M:7</td>
      <td>2</td>
      <td>1373.560</td>
      <td>8.298608e-05</td>
      <td>0</td>
      <td>False</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
<p>27 rows × 17 columns</p>
</div>


Here is a way to get the hits with reverse protein sequences. As I mentioned earlier, if you use the option --append with the perl script, you will end up with protein sequences that have ###REV### at the beginning appended for a decoy protein sequence.

```python
# Probabaly point a link to decoy.pl script?
omssa_output[' IsReverse'] = omssa_output[' Defline'].str.match('###REV###')
```


```python
omssa_output[' IsReverse'].head()
```




    0    False
    1    False
    2    False
    3    False
    4    False
    Name:  IsReverse, dtype: bool




```python
Reverse_yeast_hits_only = omssa_output[omssa_output[' IsReverse'].values]
```


```python
# All reverse hits seem to be in 0 to -2 region, which is expected.
eps=1e-32 # To take care of zeoes
fig, ax = plt.subplots()
plt.hist(np.log10(eps + Reverse_yeast_hits_only[' E-value']))
loc = ticker.MultipleLocator(base=1) # thanks again to google!
ax.xaxis.set_major_locator(loc)
ax.set_xlabel('Score in E-value')
ax.set_ylabel('Counts')
ax.set_title('Reverse hits distribution')
fig.tight_layout
```




    <bound method Figure.tight_layout of <Figure size 432x288 with 1 Axes>>




![png](output_13_1.png)

Let's remove the reverse protein sequences and focus only on the target hits.

```python
Target_yeast_hits_only = omssa_output[~omssa_output[' IsReverse'].values] # See the tilde sign
```


```python
Target_yeast_hits_only.reset_index(drop=True)
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
      <th>Spectrum number</th>
      <th>Filename/id</th>
      <th>Peptide</th>
      <th>E-value</th>
      <th>Mass</th>
      <th>gi</th>
      <th>Accession</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Defline</th>
      <th>Mods</th>
      <th>Charge</th>
      <th>Theo Mass</th>
      <th>P-value</th>
      <th>NIST score</th>
      <th>IsReverse</th>
      <th>Length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>120</td>
      <td>index=120</td>
      <td>GSIDEQHPR</td>
      <td>0.099394</td>
      <td>1037.487</td>
      <td>0</td>
      <td>BL_ORD_ID:50</td>
      <td>250</td>
      <td>258</td>
      <td>sp|P06169|PDC1_YEAST Pyruvate decarboxylase is...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1037.490</td>
      <td>3.803824e-05</td>
      <td>0</td>
      <td>False</td>
      <td>9</td>
    </tr>
    <tr>
      <th>1</th>
      <td>120</td>
      <td>index=120</td>
      <td>GSIDEQHPR</td>
      <td>0.099394</td>
      <td>1037.487</td>
      <td>0</td>
      <td>BL_ORD_ID:131</td>
      <td>250</td>
      <td>258</td>
      <td>sp|P26263|PDC6_YEAST Pyruvate decarboxylase is...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1037.490</td>
      <td>3.803824e-05</td>
      <td>0</td>
      <td>False</td>
      <td>9</td>
    </tr>
    <tr>
      <th>2</th>
      <td>339</td>
      <td>index=339</td>
      <td>TSGRPIKGDSSAGGK</td>
      <td>0.001174</td>
      <td>1416.728</td>
      <td>0</td>
      <td>BL_ORD_ID:2647</td>
      <td>176</td>
      <td>190</td>
      <td>sp|P47075|VTC4_YEAST Vacuolar transporter chap...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1416.731</td>
      <td>5.359808e-07</td>
      <td>0</td>
      <td>False</td>
      <td>15</td>
    </tr>
    <tr>
      <th>3</th>
      <td>421</td>
      <td>index=421</td>
      <td>TSGRPIKGDSSAGGK</td>
      <td>0.035522</td>
      <td>1416.729</td>
      <td>0</td>
      <td>BL_ORD_ID:2647</td>
      <td>176</td>
      <td>190</td>
      <td>sp|P47075|VTC4_YEAST Vacuolar transporter chap...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1416.731</td>
      <td>1.624232e-05</td>
      <td>0</td>
      <td>False</td>
      <td>15</td>
    </tr>
    <tr>
      <th>4</th>
      <td>506</td>
      <td>index=506</td>
      <td>NEETSGEGGEDKNEPSSK</td>
      <td>0.028961</td>
      <td>1892.786</td>
      <td>0</td>
      <td>BL_ORD_ID:4527</td>
      <td>76</td>
      <td>93</td>
      <td>sp|Q02776|TIM50_YEAST Mitochondrial import inn...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1892.789</td>
      <td>1.991815e-05</td>
      <td>0</td>
      <td>False</td>
      <td>18</td>
    </tr>
    <tr>
      <th>5</th>
      <td>514</td>
      <td>index=514</td>
      <td>ANNSQESNNATSSTSQGTR</td>
      <td>0.255101</td>
      <td>1952.833</td>
      <td>0</td>
      <td>BL_ORD_ID:1296</td>
      <td>748</td>
      <td>766</td>
      <td>sp|P39001|UME6_YEAST Transcriptional regulator...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1952.844</td>
      <td>1.608454e-04</td>
      <td>0</td>
      <td>False</td>
      <td>19</td>
    </tr>
    <tr>
      <th>6</th>
      <td>515</td>
      <td>index=515</td>
      <td>NLSDEKNDSR</td>
      <td>0.234685</td>
      <td>1176.535</td>
      <td>0</td>
      <td>BL_ORD_ID:1997</td>
      <td>363</td>
      <td>372</td>
      <td>sp|P0C2I8|YL14A_YEAST Transposon Ty1-LR4 Gag p...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1176.538</td>
      <td>8.779853e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>7</th>
      <td>515</td>
      <td>index=515</td>
      <td>NLSDEKNDSR</td>
      <td>0.234685</td>
      <td>1176.535</td>
      <td>0</td>
      <td>BL_ORD_ID:3619</td>
      <td>363</td>
      <td>372</td>
      <td>sp|P0C2J1|YP14B_YEAST Transposon Ty1-PR3 Gag-P...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1176.538</td>
      <td>8.779853e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>8</th>
      <td>515</td>
      <td>index=515</td>
      <td>NLSDEKNDSR</td>
      <td>0.234685</td>
      <td>1176.535</td>
      <td>0</td>
      <td>BL_ORD_ID:2653</td>
      <td>363</td>
      <td>372</td>
      <td>sp|P0CX71|YE11A_YEAST Transposon Ty1-ER1 Gag p...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1176.538</td>
      <td>8.779853e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>9</th>
      <td>515</td>
      <td>index=515</td>
      <td>NLSDEKNDSR</td>
      <td>0.234685</td>
      <td>1176.535</td>
      <td>0</td>
      <td>BL_ORD_ID:2713</td>
      <td>363</td>
      <td>372</td>
      <td>sp|P0C2I6|YL13B_YEAST Transposon Ty1-LR3 Gag-P...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1176.538</td>
      <td>8.779853e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
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
    </tr>
    <tr>
      <th>1670</th>
      <td>9973</td>
      <td>index=9973</td>
      <td>QSESEKLEYK</td>
      <td>0.032187</td>
      <td>1239.595</td>
      <td>0</td>
      <td>BL_ORD_ID:4272</td>
      <td>133</td>
      <td>142</td>
      <td>sp|Q06511|RRP15_YEAST Ribosomal RNA-processing...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1239.600</td>
      <td>1.431807e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>1671</th>
      <td>9974</td>
      <td>index=9974</td>
      <td>RYEAYSQQMK</td>
      <td>0.034240</td>
      <td>1302.597</td>
      <td>0</td>
      <td>BL_ORD_ID:3352</td>
      <td>776</td>
      <td>785</td>
      <td>sp|P25694|CDC48_YEAST Cell division control pr...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1302.603</td>
      <td>1.528556e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>1672</th>
      <td>9981</td>
      <td>index=9981</td>
      <td>AKDEEAAFAR</td>
      <td>0.173263</td>
      <td>1106.533</td>
      <td>0</td>
      <td>BL_ORD_ID:5971</td>
      <td>1008</td>
      <td>1017</td>
      <td>sp|P07702|LYS2_YEAST L-2-aminoadipate reductas...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1106.536</td>
      <td>6.739136e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>1673</th>
      <td>9982</td>
      <td>index=9982</td>
      <td>NSLDSSIDKK</td>
      <td>0.075037</td>
      <td>1106.533</td>
      <td>0</td>
      <td>BL_ORD_ID:5785</td>
      <td>203</td>
      <td>212</td>
      <td>sp|P09032|EI2BG_YEAST Translation initiation f...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1105.562</td>
      <td>2.918573e-05</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>1674</th>
      <td>9986</td>
      <td>index=9986</td>
      <td>NAEQAVSAEER</td>
      <td>0.512250</td>
      <td>1202.551</td>
      <td>0</td>
      <td>BL_ORD_ID:237</td>
      <td>389</td>
      <td>399</td>
      <td>sp|Q07896|NOC3_YEAST Nucleolar complex-associa...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1202.554</td>
      <td>2.276667e-04</td>
      <td>0</td>
      <td>False</td>
      <td>11</td>
    </tr>
    <tr>
      <th>1675</th>
      <td>9987</td>
      <td>index=9987</td>
      <td>KVTQMTPAPK</td>
      <td>0.812181</td>
      <td>1100.573</td>
      <td>0</td>
      <td>BL_ORD_ID:552</td>
      <td>16</td>
      <td>25</td>
      <td>sp|O14455|RL36B_YEAST 60S ribosomal protein L3...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1099.607</td>
      <td>3.103480e-04</td>
      <td>0</td>
      <td>False</td>
      <td>10</td>
    </tr>
    <tr>
      <th>1676</th>
      <td>9988</td>
      <td>index=9988</td>
      <td>KLDELENK</td>
      <td>0.006235</td>
      <td>987.523</td>
      <td>0</td>
      <td>BL_ORD_ID:6104</td>
      <td>482</td>
      <td>489</td>
      <td>sp|P31539|HS104_YEAST Heat shock protein 104 O...</td>
      <td>NaN</td>
      <td>2</td>
      <td>987.525</td>
      <td>2.068072e-06</td>
      <td>0</td>
      <td>False</td>
      <td>8</td>
    </tr>
    <tr>
      <th>1677</th>
      <td>9988</td>
      <td>index=9988</td>
      <td>KVEELENK</td>
      <td>0.138069</td>
      <td>987.523</td>
      <td>0</td>
      <td>BL_ORD_ID:2693</td>
      <td>109</td>
      <td>116</td>
      <td>sp|Q12172|YRR1_YEAST Zinc finger transcription...</td>
      <td>NaN</td>
      <td>2</td>
      <td>987.525</td>
      <td>4.579411e-05</td>
      <td>0</td>
      <td>False</td>
      <td>8</td>
    </tr>
    <tr>
      <th>1678</th>
      <td>9995</td>
      <td>index=9995</td>
      <td>LGSSSQSFK</td>
      <td>0.025725</td>
      <td>939.463</td>
      <td>0</td>
      <td>BL_ORD_ID:5255</td>
      <td>222</td>
      <td>230</td>
      <td>sp|P32469|DPH5_YEAST Diphthine methyl ester sy...</td>
      <td>NaN</td>
      <td>2</td>
      <td>939.466</td>
      <td>1.216335e-05</td>
      <td>0</td>
      <td>False</td>
      <td>9</td>
    </tr>
    <tr>
      <th>1679</th>
      <td>9996</td>
      <td>index=9996</td>
      <td>VTSSTSSSSSNGSVDVR</td>
      <td>0.000006</td>
      <td>1655.753</td>
      <td>0</td>
      <td>BL_ORD_ID:644</td>
      <td>94</td>
      <td>110</td>
      <td>sp|P40335|PEP8_YEAST Carboxypeptidase Y-defici...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1655.759</td>
      <td>3.424696e-09</td>
      <td>0</td>
      <td>False</td>
      <td>17</td>
    </tr>
  </tbody>
</table>
<p>1680 rows × 17 columns</p>
</div>


This is the score distribution, but using only target hits.

```python
# All reverse hits seem to be in 0 to -2 region, which is expected.
eps=1e-32 # To take care of zeoes
fig, ax = plt.subplots()
plt.hist(np.log10(eps + Target_yeast_hits_only[' E-value']), bins = 100)
loc = ticker.MultipleLocator(base=2) # thanks again to google!
ax.xaxis.set_major_locator(loc)
ax.set_xlabel('Score in E-value')
ax.set_ylabel('Counts')
ax.set_title('Target hits distribution')
fig.tight_layout
```




    <bound method Figure.tight_layout of <Figure size 432x288 with 1 Axes>>




![png](output_16_1.png)



```python
Reverse_yeast_hits_only.head()
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
      <th>Spectrum number</th>
      <th>Filename/id</th>
      <th>Peptide</th>
      <th>E-value</th>
      <th>Mass</th>
      <th>gi</th>
      <th>Accession</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Defline</th>
      <th>Mods</th>
      <th>Charge</th>
      <th>Theo Mass</th>
      <th>P-value</th>
      <th>NIST score</th>
      <th>IsReverse</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>75</th>
      <td>872</td>
      <td>index=872</td>
      <td>QETASNKPLKLYSCITR</td>
      <td>0.456971</td>
      <td>2007.920</td>
      <td>0</td>
      <td>BL_ORD_ID:9855</td>
      <td>835</td>
      <td>851</td>
      <td>###REV###sp|Q08387|DNLI4_YEAST Reverse sequenc...</td>
      <td>NaN</td>
      <td>3</td>
      <td>2008.042</td>
      <td>0.000310</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>117</th>
      <td>1103</td>
      <td>index=1103</td>
      <td>RTLEPTSLGGLIEVLR</td>
      <td>0.574082</td>
      <td>1753.787</td>
      <td>0</td>
      <td>BL_ORD_ID:11219</td>
      <td>659</td>
      <td>674</td>
      <td>###REV###sp|P38850|RT107_YEAST Reverse sequenc...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1753.010</td>
      <td>0.000363</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>131</th>
      <td>1179</td>
      <td>index=1179</td>
      <td>KASLLILDDHSDDNK</td>
      <td>0.381089</td>
      <td>1682.339</td>
      <td>0</td>
      <td>BL_ORD_ID:10482</td>
      <td>630</td>
      <td>644</td>
      <td>###REV###sp|Q06673|ECM30_YEAST Reverse sequenc...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1682.848</td>
      <td>0.000249</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>167</th>
      <td>1375</td>
      <td>index=1375</td>
      <td>KASLLILDDHSDDNK</td>
      <td>0.378321</td>
      <td>1683.785</td>
      <td>0</td>
      <td>BL_ORD_ID:10482</td>
      <td>630</td>
      <td>644</td>
      <td>###REV###sp|Q06673|ECM30_YEAST Reverse sequenc...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1682.848</td>
      <td>0.000248</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>595</th>
      <td>3081</td>
      <td>index=3081</td>
      <td>DAAQVAEEVDDER</td>
      <td>0.733552</td>
      <td>1443.671</td>
      <td>0</td>
      <td>BL_ORD_ID:11559</td>
      <td>428</td>
      <td>440</td>
      <td>###REV###sp|P46367|ALDH4_YEAST Reverse sequenc...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1445.628</td>
      <td>0.000367</td>
      <td>0</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>




```python
Target_yeast_hits_only.head()
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
      <th>Spectrum number</th>
      <th>Filename/id</th>
      <th>Peptide</th>
      <th>E-value</th>
      <th>Mass</th>
      <th>gi</th>
      <th>Accession</th>
      <th>Start</th>
      <th>Stop</th>
      <th>Defline</th>
      <th>Mods</th>
      <th>Charge</th>
      <th>Theo Mass</th>
      <th>P-value</th>
      <th>NIST score</th>
      <th>IsReverse</th>
      <th>Length</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>120</td>
      <td>index=120</td>
      <td>GSIDEQHPR</td>
      <td>0.099394</td>
      <td>1037.487</td>
      <td>0</td>
      <td>BL_ORD_ID:50</td>
      <td>250</td>
      <td>258</td>
      <td>sp|P06169|PDC1_YEAST Pyruvate decarboxylase is...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1037.490</td>
      <td>3.803824e-05</td>
      <td>0</td>
      <td>False</td>
      <td>9</td>
    </tr>
    <tr>
      <th>1</th>
      <td>120</td>
      <td>index=120</td>
      <td>GSIDEQHPR</td>
      <td>0.099394</td>
      <td>1037.487</td>
      <td>0</td>
      <td>BL_ORD_ID:131</td>
      <td>250</td>
      <td>258</td>
      <td>sp|P26263|PDC6_YEAST Pyruvate decarboxylase is...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1037.490</td>
      <td>3.803824e-05</td>
      <td>0</td>
      <td>False</td>
      <td>9</td>
    </tr>
    <tr>
      <th>2</th>
      <td>339</td>
      <td>index=339</td>
      <td>TSGRPIKGDSSAGGK</td>
      <td>0.001174</td>
      <td>1416.728</td>
      <td>0</td>
      <td>BL_ORD_ID:2647</td>
      <td>176</td>
      <td>190</td>
      <td>sp|P47075|VTC4_YEAST Vacuolar transporter chap...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1416.731</td>
      <td>5.359808e-07</td>
      <td>0</td>
      <td>False</td>
      <td>15</td>
    </tr>
    <tr>
      <th>3</th>
      <td>421</td>
      <td>index=421</td>
      <td>TSGRPIKGDSSAGGK</td>
      <td>0.035522</td>
      <td>1416.729</td>
      <td>0</td>
      <td>BL_ORD_ID:2647</td>
      <td>176</td>
      <td>190</td>
      <td>sp|P47075|VTC4_YEAST Vacuolar transporter chap...</td>
      <td>NaN</td>
      <td>2</td>
      <td>1416.731</td>
      <td>1.624232e-05</td>
      <td>0</td>
      <td>False</td>
      <td>15</td>
    </tr>
    <tr>
      <th>4</th>
      <td>506</td>
      <td>index=506</td>
      <td>NEETSGEGGEDKNEPSSK</td>
      <td>0.028961</td>
      <td>1892.786</td>
      <td>0</td>
      <td>BL_ORD_ID:4527</td>
      <td>76</td>
      <td>93</td>
      <td>sp|Q02776|TIM50_YEAST Mitochondrial import inn...</td>
      <td>NaN</td>
      <td>3</td>
      <td>1892.789</td>
      <td>1.991815e-05</td>
      <td>0</td>
      <td>False</td>
      <td>18</td>
    </tr>
  </tbody>
</table>
</div>




```python
Target_yeast_hits_only.shape
```




    (1680, 17)




```python
Target_yeast_hits_only[' Length'].head()
```




    0     9
    1     9
    2    15
    3    15
    4    18
    Name:  Length, dtype: int64


Earlier I mentioned that the program scores a hit if the observed peptide mass matches the experimental. We consider a small offset, as the mass cannot be exactly the same. This offset is generally in ppm (parts per million), and in most of the new high-res instruments, this value is approximately 20ppm.  

20ppm at 5000Da is 0.1. You can expect 0.1 difference at 5K peptide mass. In the below plot, you can see that the distribution is very tight around 0.


```python
fig, ax = plt.subplots()
plt.hist(Target_yeast_hits_only[' Mass'] - Target_yeast_hits_only[' Theo Mass'], bins = 5000)
loc = ticker.MultipleLocator(base=0.01) # thanks again to google!
ax.xaxis.set_major_locator(loc)
ax.set_xlabel('delta Precursor Mass [Da]')
ax.set_ylabel('Counts')
ax.set_xlim([-0.02, 0.02])
ax.set_title("Difference of observed and expected peptide mass")
fig.tight_layout
```




    <bound method Figure.tight_layout of <Figure size 432x288 with 1 Axes>>




![png](output_21_1.png)

Let's take a look at the length of the peptides present in the sample (i.e., MGF file)

```python
Target_yeast_hits_only[' Charge'].astype('category')
fig, ax = plt.subplots()
plt.hist(Target_yeast_hits_only[' Length'], bins = 100)
loc = ticker.MultipleLocator(base=5) # thanks again to google!
ax.xaxis.set_major_locator(loc)
ax.set_xlabel('Length')
ax.set_ylabel('Counts')
fig.tight_layout
```




    <bound method Figure.tight_layout of <Figure size 432x288 with 1 Axes>>




![png](output_22_1.png)

Let's take a look at the charge of the peptides present in the sample (i.e., MGF file)

```python
fig, ax = plt.subplots()
plt.hist(Target_yeast_hits_only[' Charge'])
loc = ticker.MultipleLocator(base=1) # this locator puts ticks at regular intervals
ax.xaxis.set_major_locator(loc)
ax.set_xlabel('Charge')
ax.set_ylabel('Counts')
fig.tight_layout
```




    <bound method Figure.tight_layout of <Figure size 432x288 with 1 Axes>>




![png](output_23_1.png)



```python

```
