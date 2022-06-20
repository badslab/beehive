SHELL=/bin/bash

.PHONY:
.ONESHELL:
fix_templates:
	for bdir in ./bokeh/*; do
		tdir="$$bdir/templates"
		echo $$bdir $$tdir
		echo "Make templates symlink"
		( cd $$bdir ; ln -sf ../../templates . )
	done

.PHONY:
.ONEHSELL:
fix_bokeh_static_js:
	JSPATH=$$(python -c 'import bokeh; print(bokeh.util.paths.bokehjsdir())')
	echo $$JSPATH
	mkdir -p static/bokeh/
	rsync  -arv $$JSPATH/ static/bokeh/
	chmod -R a+rX static/bokeh


.ONESHELL:
serve_cbd2: fix_templates fix_bokeh_static_js
	while true; do 
		echo "(re)starting)"
		bokeh serve --use-xheaders \
			--allow-websocket-origin=data.bdslab.org \
			 --port 5009 bokeh/gene_expression/
		sleep 0.5
	done

.ONEHSELL:
serve_moam: fix_templates
	pipenv run bokeh serve --port 5010 bokeh/gene_expression/  bokeh/diffexp/

.ONEHSELL:
serve_sjo: fix_templates
	bokeh serve --port 5010 bokeh/gene_expression/  bokeh/diffexp/

.ONEHSELL:
serve_sjo_dev: fix_templates
	bokeh serve --dev --port 5011 bokeh/gene_expression/

.PHONY:
.ONESHELL:
static_website:
	cd static
	make html
	git add content/*.md output/*.html  output/category/*.html output/author/*.html output/tag/*.html
	git commit -m 'rebuild static website' content/*.md output/*.html  output/category/*.html output/author/*.html output/tag/*.html


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
