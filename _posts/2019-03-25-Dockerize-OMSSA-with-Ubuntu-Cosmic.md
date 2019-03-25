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


Writing workflows and pipelines in proteomics informatics is still emerging. Proteome Discoverer is one such software where you can embed different software and edit different sections of the pipeline/workflow, depending on the needs. Here, I will NOT focus on any such workflows, but plan to use docker containers to show that such pipelines can be possible, but might require a significant effort.

With so many variables at different parts of the pipeline/workflow, reproducible research might be an issue, if proper caution is not taken care of. Here, I will introduce briefly to the idea of using dockerfile with a proteomics search engine. I did not “google”, but I have not seen any work using Docker containers in proteomics pipelines. To keep things simple, I did an OMSSA search on a sample data file within a container. It is possible to extend this analyses to have similar Python/R scripts to get protein lists and further extend the analyses.

For the current blog post, I use linux version of OMSSA that can be found [here](ftp://ftp.ncbi.nih.gov/pub/lewisg/omssa/CURRENT/). I use ubuntu:cosmic within the docker. For makeblastdb, we need ncbi-blast+ found [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

The end goal of this blog is a docker container that you can use for your proteomics research with (1) your input mgf file and (2) your fasta file of interest. This dockerfile will help to set-up your OMSSA search within few minutes.

However if you use this container, also NOTE to edit the omssa_run.sh file with the arguments that work best for the OMSSA database search (for example, add phosphorylation as variable modification if you are searching a phosphoproteomics dataset).

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

I do not want to go through the details of the dockerfile, but all I did was to get the tar.gz of OMSSA and ncbi-blast+, copy sample data sets to the container and run makeblastdb followed by an OMSSA search.
