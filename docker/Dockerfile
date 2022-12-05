# syntax=docker/dockerfile:1
FROM python:3.8

WORKDIR beehive/
COPY . .

RUN set -x \
    # create nginx user/group first, to be consistent throughout docker variants (borrowed from nginx docker image)
    && addgroup --system  nginx \
    && adduser --system --disabled-login --ingroup nginx --no-create-home --home /nonexistent --gecos "nginx user" --shell /bin/false nginx \
    && apt-get update \
    && apt-get install --no-install-recommends --no-install-suggests -y zile gnupg1 ca-certificates build-essential nginx \
    && pip install --upgrade pip polars pyarrow colorcet bokeh==2.4.3 matplotlib scanpy pyyaml pandas scikit-learn \
    && pip install --upgrade scipy statsmodels typer pelican markdown pymed \
    && pip install -e beehive \
    && chmod a+x docker/entrypoint.sh \
    o    && make fix_templates \
    && make fix_bokeh_static_js

ENV BEEHIVE_BASEDIR=/beehive/data

# nginx reverse proxy
EXPOSE 5010

ENTRYPOINT ["./docker/entrypoint.sh"]