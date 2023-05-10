#!/bin/bash

set -v

echo "Starting Docker Entrypoint"

# start nginx to serve static pages & as a reverse proxy
echo "Start nginx"
nginx -c /beehive/docker/nginx.conf


echo "Go in infinite bokeh loop"
while true; do
    echo "Starting bokeh"
    bokeh serve \
    --allow-websocket-origin=${VISUALIZATION_WEBSOCKET_ORIGIN} \
    --port 5009 \
    bokeh/gene_expression/ \
    bokeh/volcano_plot/ \
    bokeh/scatter_expression/ \
    bokeh/hexbin_expression/ \
    bokeh/quadrant_plot/ \
    bokeh/jitter_expression/ \
    bokeh/gsea_view/
    echo "Crashed? That's not good..."
    echo "-------------------------------------"
    sleep 1
done