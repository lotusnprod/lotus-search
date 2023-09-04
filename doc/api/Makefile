.PHONY: doc lint serve
doc: lint
	mkdir -p dist
	docker run --rm -e REDOCLY_TELEMETRY=off -v $(CURDIR):/spec redocly/cli bundle -o dist/lotus_api.yaml lotus_api.yaml
	cp rapidoc/index.html dist/index.html
	echo "You can see the documentation by running 'make serve'"
lint:
	docker run --rm -e REDOCLY_TELEMETRY=off -v $(CURDIR):/spec redocly/cli lint
serve: doc
	docker run --rm --name lotus_api_doc -p 8077:80 -v $(CURDIR)/dist:/usr/share/nginx/html:ro nginx:alpine
