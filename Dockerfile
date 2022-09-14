# syntax=docker/dockerfile:1
FROM python:3.8
WORKDIR beehive/
COPY . .
ENV BEEHIVE_BASEDIR=/beehive
COPY ./requirements.txt ./requirements.txt
RUN pip install -r ./requirements.txt --default-timeout=1000
# COPY ./data ./data
# RUN mkdir ./data/h5ad
# RUN ln -s ./data/h5ad_store/*.yaml ./data/h5ad/
COPY ./Makefile ./Makefile
RUN make fix_templates