SHELL=/bin/bash

.PHONY:
.SILENT:
check_deployment:
	echo "Pre deployment check in $$(hostname):$$PWD"
	[ -d ./data/h5ad ] || (echo "No ./data/h5ad folder - please create " ; false)
	(ls ./data/h5ad/*yaml > /dev/null ) ||  (echo "No yaml files in data/h5ad folder - please populate" ; false)
	for file in ./data/h5ad/*yaml; do \
		bn=$$(basename $$file .yaml); \
		[ -f ./data/h5ad/$${bn}.X.prq ] || (echo "Can't find $${bn}.X.prq"; false); \
		[ -f ./data/h5ad/$${bn}.obs.prq ] || (echo "Can't find $${bn}.obs.prq"; false); \
		[ -f ./data/h5ad/$${bn}.var.prq ] || (echo "Can't find $${bn}.var.prq"; false); \
	done



.PHONY:
fix_templates:
	for bdir in ./bokeh/*; do \
		tdir="$$bdir/templates"; \
		echo $$bdir $$tdir; \
		echo "Make templates symlink"; \
		( cd $$bdir ; ln -sf ../../templates . ); \
	done

.PHONY:
.ONEHSELL:
fix_bokeh_static_js:
	JSPATH=$$(python -c 'import bokeh; print(bokeh.util.paths.bokehjsdir())')
	echo "Getting bokeh static data from: $$JSPATH"
	mkdir -p static/bokeh/
	rsync  -arv $$JSPATH/ static/bokeh/
	chmod -R a+rX static/bokeh
	mkdir -p static/bokeh/
	rsync  -arv $$JSPATH/ static/bokeh/
	chmod -R a+rX static/bokeh


.ONESHELL:
serve_cbd2: fix_templates fix_bokeh_static_js check_deployment
	export BEEHIVE_BASEDIR=/media/gbw_cbdbads_alzmap/bdslab_visualization_public/beehive
	while true; do \
		echo "(re)starting)"		
		bokeh serve --use-xheaders \
			--allow-websocket-origin=data.bdslab.org \
			 --port 5009 \
			bokeh/gene_expression/ \
			bokeh/volcano_plot/ \
			bokeh/scatter_expression/
		sleep 0.5
	done

.ONESHELL:
serve_cbd2_private: fix_templates fix_bokeh_static_js check_deployment
	export BEEHIVE_BASEDIR=/media/gbw_cbdbads_alzmap/bdslab_visualization_private/beehive
	while true; do 
		echo "(re)starting)" 
		bokeh serve --use-xheaders \
			--allow-websocket-origin=muna.bdslab.org \
			 --port 5008 \
			bokeh/gene_expression/ \
			bokeh/volcano_plot/ \
			bokeh/scatter_expression/
		sleep 0.5
	done

.ONEHSELL:
serve_dev: fix_templates
	bokeh serve --dev --port 5010 bokeh/scatter_expression/ # bokeh/diffexp/
	#	bokeh serve --dev --port 5010 bokeh/gene_expression/ # bokeh/diffexp/


serve_dev_raghid: fix_templates
	pipenv run bokeh serve --dev --port 5009 --allow-websocket-origin=* bokeh/gene_expression/

serve_dev_raghid_auth: fix_templates
	pipenv run bokeh serve --dev --port 5009 --allow-websocket-origin=* --auth-module=beehive/beehive/auth.py bokeh/gene_expression/

serve_dev_raghid_scatter_local: fix_templates
	pipenv run bokeh serve --dev --port 5009 --allow-websocket-origin=* --auth-module=beehive/beehive/auth.py --ssl-certfile beehive/cert.pem --ssl-keyfile beehive/key.pem bokeh/scatter_expression/

serve_dev_raghid_scatter_server: fix_templates
	pipenv run bokeh serve --dev --port 8009 --allow-websocket-origin=* --auth-module=beehive/beehive/auth.py --ssl-certfile beehive/cert.pem --ssl-keyfile beehive/key.pem bokeh/scatter_expression/

serve_dev_raghid2: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/scatter_expression/

serve_dev_raghid3: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/hexbin_expression/

serve_dev_raghid4: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/volcano_plot/



.PHONY:
.ONESHELL:
rebuild_static_private_website:
	cd static_private
	make html
	git add content/*.md output/*.html  output/category/*.html output/author/*.html output/tag/*.html
	git commit -m 'rebuild private static website' content/*.md output/*.html  output/category/*.html output/author/*.html output/tag/*.html


.PHONY:
.ONESHELL:
rebuild_static_public_website:
	cd static_public
	make html
	git add content/*.md output/*.html  output/category/*.html output/author/*.html output/tag/*.html
	git commit -m 'rebuild static website' content/*.md output/*.html  output/category/*.html output/author/*.html output/tag/*.html

sync_data_to_moamx:
	rsync -arv data/h5ad/*prq moamx:/data/project/mark/beehive/data/h5ad/


# .SILENT:
# .ONESHELL:
# dev:
# 	export PATH=/opt/python/miniconda-sep-2021/envs/sparrowhawk/bin/:$$PATH;
# 	LATEST=$$(ls -dt bokeh/*/main.py | head -1);
# 	LATEST=$$(readlink -f $$LATEST);
# 	LATEST=$$(dirname $$LATEST);
# 	echo "Using bokeh app: $$LATEST"
# 	while true; do
# 		echo "(re-)starting $$LATEST"
# 		pipenv run bokeh serve --dev --port 5008 $$LATEST;
# 	done
