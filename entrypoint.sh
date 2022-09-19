#!/bin/bash
bokeh serve --allow-websocket-origin=${VISUALIZATION_WEBSOCKET_ORIGIN} \
			 --port ${PORT_FOR_VISUALIZATION} \
			bokeh/gene_expression/ \
			bokeh/volcano_plot/ \
			bokeh/scatter_expression/ \
            bokeh/hexbin_expression/ \
            bokeh/quadrant_plot/ 