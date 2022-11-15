# syntax=docker/dockerfile:1
FROM python:3.8
WORKDIR beehive/
COPY . .
RUN apt-get update
RUN apt-get install -y build-essential
RUN  pip install --upgrade pip polars pyarrow colorcet bokeh==2.4.3 matplotlib scanpy pyyaml pandas scikit-learn scipy statsmodels typer pelican markdown pymed
RUN pip install -e beehive
COPY ./Makefile ./Makefile
ENV BEEHIVE_BASEDIR=/beehive/data
EXPOSE ${PORT_FOR_VISUALIZATION}
RUN make fix_templates
RUN ["chmod", "+x", "./entrypoint.sh"]
ENTRYPOINT ["./entrypoint.sh"]