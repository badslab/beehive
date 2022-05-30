SHELL=/bin/bash

.ONEHSELL:
serve_moam:
	pipenv run bokeh serve --port 5009 bokeh/gene_expression/


.SILENT:
.ONESHELL:
dev:
	export PATH=/opt/python/miniconda-sep-2021/envs/sparrowhawk/bin/:$$PATH;
	LATEST=$$(ls -dt bokeh/*/main.py | head -1);
	LATEST=$$(readlink -f $$LATEST);
	LATEST=$$(dirname $$LATEST);
	echo "Using bokeh app: $$LATEST"
	while true; do
		echo "(re-)starting $$LATEST"
		pipenv run bokeh serve --dev --port 5008 $$LATEST;
	done
