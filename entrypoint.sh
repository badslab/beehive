#!/bin/bash
bokeh serve --allow-websocket-origin=* \
			 --port 5009 \
			bokeh/gene_expression/ \
			bokeh/volcano_plot/ \
			bokeh/scatter_expression/ \
            bokeh/hexbin_expression/ \
            bokeh/quadrant_plot/ 