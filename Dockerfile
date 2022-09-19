# syntax=docker/dockerfile:1
FROM python:3.8
RUN apt-get -y update
RUN apt-get -y install git
WORKDIR beehive/
COPY . .
ENV BEEHIVE_BASEDIR=/beehive
COPY ./requirements.txt ./requirements.txt
RUN pip install -r ./requirements.txt --default-timeout=1000
COPY ./Makefile ./Makefile
EXPOSE ${PORT_FOR_VISUALIZATION}
ENV PORT_FOR_VISUALIZATION=${PORT_FOR_VISUALIZATION}
ENV VISUALIZATION_WEBSOCKET_ORIGIN=${VISUALIZATION_WEBSOCKET_ORIGIN}
RUN make fix_templates
RUN ["chmod", "+x", "./entrypoint.sh"]
ENTRYPOINT ["./entrypoint.sh"]