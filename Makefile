

.ONEHSELL:
serve_moam:
	export PATH=/opt/python/miniconda-sep-2021/envs/sparrowhawk/bin/:$PATH
	bokeh serve --port 5009 bokeh/gene_expression/
