SHELL=/bin/bash

fix_templates:
	for bdir in ./bokeh/*; do
		tdir="$$bdir/templates"
		echo $$bdir $$tdir
		echo "Make templates symlink"
		( cd $$bdir ; ln -sf ../../templates . )
	done

.ONESHELL:
serve_cbd2: fix_templates
	while true; do 
		echo "(re)starting)"
		bokeh serve --use-xheaders \
			--allow-websocket-origin=data.bdslab.org \
			 --port 5009 bokeh/gene_expression/	
		sleep 0.5
	done

.ONEHSELL:
serve_moam: fix_templates
	pipenv run bokeh serve --port 5009 bokeh/gene_expression/  bokeh/diffexp/

.ONEHSELL:
serve_sjo: fix_templates
	bokeh serve --port 5009 bokeh/gene_expression/  bokeh/diffexp/

.ONEHSELL:
serve_sjo_dev: fix_templates
	bokeh serve --dev --port 5009 bokeh/gene_expression/

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
