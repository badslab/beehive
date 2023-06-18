# syntax=docker/dockerfile:1
FROM python:3.8

WORKDIR /beehive/
COPY . .

RUN set -x \
    # create nginx user/group first, to be consistent throughout docker variants (borrowed from nginx docker image)
    && addgroup --system  nginx \
    && adduser --system --disabled-login --ingroup nginx --no-create-home --home /nonexistent --gecos "nginx user" --shell /bin/false nginx \
    && apt-get update \
    && apt-get install -y --no-install-recommends --no-install-suggests zile gnupg1 less \
    && apt-get install -y --no-install-recommends --no-install-suggests ca-certificates build-essential nginx \
    && pip install --upgrade pip wheel \
    && pip install --upgrade pip polars pyarrow colorcet bokeh==2.4.3 matplotlib scanpy pyyaml h5py \
    && pip install --upgrade markdown pandas scikit-learn scipy statsmodels typer pelican pymed \
    && pip install -e beehive \
    && chmod a+x docker/entrypoint.sh \
    && make fix_templates \
    && make fix_bokeh_static_js

ENV BEEHIVE_BASEDIR=/beehive/data

# nginx reverse proxy
EXPOSE 5010

ENTRYPOINT ["./docker/entrypoint.sh"]