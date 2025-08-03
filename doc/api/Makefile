.PHONY: image doc lint serve ts
image:
	docker build -t openapi_tools .
doc: image
	mkdir -p dist
	docker run --rm -e REDOCLY_TELEMETRY=off -v $(CURDIR):/spec redocly/cli bundle -o dist/lotus_api.yaml lotus_api.yaml
	cp rapidoc/index.html dist/index.html
	echo "You can see the documentation by running 'make serve'"
lint: doc
	docker run --rm -e REDOCLY_TELEMETRY=off -v $(CURDIR):/spec redocly/cli lint
serve: image doc
	docker run --rm --name lotus_api_doc -p 8077:80 -v $(CURDIR)/dist:/usr/share/nginx/html:ro nginx:alpine
ts: image doc lint
	docker run --rm -v $(CURDIR)/dist:/spec -v $(CURDIR)/../../frontend/src/interfaces:/ts openapi_tools openapi-typescript -o /ts/lotus_api.ts /spec/lotus_api.yaml
