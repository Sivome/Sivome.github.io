---
title: "DOCKERize OMSSA with UBUNTU Cosmic"
date: '2019-03-25'
layout: post
categories: Proteomics
---


If you’re working in mass-spectrometry based proteomics informatics, you already know the effort it takes to systematically analyze the mass-spec data to generate processed data. For example,
choosing tools to convert raw data to another format,
choosing a database search fasta file for the sample of interest,
choosing sample-specific parameters (e.g., search for phosphorylation on serine/threonine/tyrosine in phospho-sample),
choosing proteomics software depending on the end goals (Identification, Quantitation).


Writing workflows and pipelines in proteomics informatics is still emerging. Proteome Discoverer is one such software where you can embed different software and edit different sections of the pipeline/workflow, depending on the needs. Here, I will NOT focus on any such workflows, but plan to use docker containers to show how docker can be used for proteomics bioinformatics to embed different tools for reproducible research.

For the current blog post, I use linux version of OMSSA that can be found [here](ftp://ftp.ncbi.nih.gov/pub/lewisg/omssa/CURRENT/). I use ubuntu:cosmic within the docker. For makeblastdb, I used ncbi-blast+ that can be found [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

The end goal of this blog is a docker container that you can use for your proteomics database search with OMSSA. Once you have your (1) your input mgf file and (2) your fasta file of interest, then you can easily tweak the below Dockerfile.

NOTE to edit the omssa_run.sh file if you want to adjust the search parameters (for example, add phosphorylation as variable modification if you are searching a phosphoproteomics dataset).

Here is the [github link](https://github.com/viswam78/dockerizeOMSSA) that has the same information and the [dockerhub link](https://cloud.docker.com/u/sridharabio/repository/docker/sridharabio/dockerize123omssa) to pull the image directly.

Here is how the dockerfile looks like:


```console
FROM ubuntu:cosmic

RUN apt-get update
RUN apt-get install -y wget

RUN wget ftp://ftp.ncbi.nih.gov/pub/lewisg/omssa/CURRENT/omssa-linux.tar.gz
RUN tar -xvf omssa-linux.tar.gz
RUN rm -rf omssa-linux.tar.gz

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz
RUN tar -xvf ncbi-blast-2.8.1+-x64-linux.tar.gz
RUN rm -rf ncbi-blast-2.8.1+-x64-linux.tar.gz

ENV PATH="${PATH}:~/../omssa-2.1.9.linux"
ENV PATH="${PATH}:~/../ncbi-blast-2.8.1+/bin"

RUN export PATH

# copy local sample files (mgf, fasta, along with mods.xml and usermods.xml)
RUN mkdir sample_data
RUN cp ../omssa-2.1.9.linux/mods.xml /sample_data
RUN cp ../omssa-2.1.9.linux/usermods.xml /sample_data
COPY S288c_RefSeq.fasta /sample_data
COPY yeast_subset.mgf /sample_data
COPY makeblastdbRUN.sh /sample_data
COPY omssa_run.sh /sample_data
```

I do not want to go through the details of the dockerfile, but all I did was [1] to get the OMSSA and ncbi-blast+ downloads, [2] to copy sample data set, along with the fasta file to the container and [3] to run makeblastdb followed by an OMSSA search.

Additional resources:
[Galaxy Proteomics](https://github.com/galaxyproteomics/docker-galaxyp)
