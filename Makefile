SHELL=/bin/bash

#BEEHIVE_BASEDIR=/Users/u0089478/data/beehive/beehive_data_intern
#BEEHIVE_BASEDIR=/Users/raghidbsat/Extras_No_Icloud_save/VIB/github_bds_lab/beehive_data_intern
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
	cp -r $$JSPATH/* static/
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
	export BEEHIVE_BASEDIR=/data/project/mark/beehive
	while true; do
		echo "(re)starting)"
		# bokeh serve --dev --port 5010 bokeh/scatter_expression/ # bokeh/diffexp/
		bokeh serve --dev --port 5010 bokeh/gene_expression/ # bokeh/diffexp/
	done

#### Docker

docker_build:
	docker buildx \
		build --platform linux/amd64 \
		-f docker/Dockerfile \
		-t bdslab/beehive \
		.

docker_build_mac:
	docker  build  -f docker/Dockerfile -t bdslab/beehive:mac_m1  .


docker_push:
	docker push bdslab/beehive:latest


docker_logs:
	export iid=$$(docker container ls -q  | head -1) ; \
	echo "checking logs of $$iid" ; \
	docker logs --follow $$iid


docker_shell:
	export iid=$$(docker container ls -q  | head -1) ; \
	echo "connecting to $$iid" ; \
	docker exec -it $$iid bash


docker_kill_mac:
	export iid=$$(docker container ls -q  | head -1) ; \
	echo "kill docker $$iid" ; \
	docker container stop -t 0 $$iid ; \
	docker container rm $$iid


.ONESHELL:
docker_run_mac:
	set -v ; \
	docker run --restart unless-stopped -d -p 5010:5010 \
		-e VISUALIZATION_WEBSOCKET_ORIGIN="*" \
		--volume ${BEEHIVE_BASEDIR}:/beehive/data \
		bdslab/beehive:mac_m1

.ONESHELL:
serve_mark_all:
	export BEEHIVE_BASEDIR=${BEEHIVE_BASEDIR} ; \
	bokeh serve \
		--port 5009 --allow-websocket-origin=*  bokeh/*


.ONESHELL:
serve_mark_gene_expression:
	export BEEHIVE_DEBUG=1 ; \
	export BEEHIVE_BASEDIR=${BEEHIVE_BASEDIR} ; \
	bokeh serve \
		--dev --port 15009 --allow-websocket-origin=* \
		bokeh/gene_expression/


.ONESHELL:
serve_mark_scatter:
	export BEEHIVE_DEBUG=1 ; \
	export BEEHIVE_BASEDIR=${BEEHIVE_BASEDIR} ; \
	while true; do bokeh serve \
		--dev --port 15009 --allow-websocket-origin=* \
		bokeh/scatter_expression/ ; sleep 1 ; done


.ONESHELL:
serve_mark_gsea:
	export BEEHIVE_DEBUG=1 ; \
	export BEEHIVE_BASEDIR=${BEEHIVE_BASEDIR} ; \
	while true; do bokeh serve \
		--dev --port 15009 --allow-websocket-origin=* \
		bokeh/gsea_view/ ; sleep 1 ; done


.ONESHELL:
serve_mark_gsea:
	export BEEHIVE_DEBUG=1 ; \
	export BEEHIVE_BASEDIR=${BEEHIVE_BASEDIR} ; \
	while true; do bokeh serve \
			--dev --port 15009 --allow-websocket-origin=* \
			bokeh/gsea_view/; done


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


##
## Targets from Raghid
##


####raghid####
#using pipenv to install requirements.txt as virtual env
serve_dev_raghid_auth: fix_templates
	pipenv run bokeh serve --dev --port 5009 --allow-websocket-origin=* --auth-module=beehive/beehive/auth.py bokeh/gene_expression/

serve_dev_raghid_auth_ssl: fix_templates
	pipenv run bokeh serve --dev --port 5009 --allow-websocket-origin=* --auth-module=beehive/beehive/auth.py --ssl-certfile beehive/cert.pem --ssl-keyfile beehive/key.pem bokeh/scatter_expression/

serve_dev_raghid_server_auth: fix_templates
	pipenv run bokeh serve --dev --port 8009 --allow-websocket-origin=* --auth-module=beehive/beehive/auth.py --ssl-certfile beehive/cert.pem --ssl-keyfile beehive/key.pem bokeh/scatter_expression/

serve_dev_raghid_gene: fix_templates
	pipenv run bokeh serve --dev --port 5009 --allow-websocket-origin=* bokeh/gene_expression/

serve_dev_raghid_scatter: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/scatter_expression/

serve_dev_raghid_hexbin: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/hexbin_expression/

serve_dev_raghid_volcano: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/volcano_plot/

serve_dev_raghid_quadrant: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/quadrant_plot/

serve_dev_raghid_abundance: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/mean_abundance/

serve_dev_raghid_all: fix_templates
	pipenv run bokeh serve --port 5009 bokeh/quadrant_plot/ \
									bokeh/volcano_plot/ \
									bokeh/scatter_expression/ \
									bokeh/gene_expression/ \
									bokeh/hexbin_expression/ \
									bokeh/jitter_expression/

serve_dev_raghid_jitter: fix_templates
	pipenv run bokeh serve --dev --port 5009 bokeh/jitter_expression/

