# syntax=docker/dockerfile:1
FROM python:3.8
WORKDIR beehive/
COPY . .
ENV BEEHIVE_BASEDIR=/beehive
COPY ./requirements.txt ./requirements.txt
RUN pip install -r ./requirements.txt --default-timeout=1000
COPY ./Makefile ./Makefile
EXPOSE 5009
RUN make fix_templates