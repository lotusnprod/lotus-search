.PHONY: all start stop test tests clean fix clean_data clean_db bench

# Standard targets
all: start

clean: stop clean_data clean_db

test: tests

# Your existing targets
start:
	docker compose up --build -d

stop:
	docker compose down

tests: start
	docker compose run -e PYTHONDONTWRITEBYTECODE=1 backend uv run pytest --cov-config=.coveragerc --cov=.

fix: start
	docker compose run backend uv run black .
	docker compose run backend uv run isort .

clean_data:
	docker compose run backend rm -rf /app/data/*

clean_db:
	docker compose run backend rm -rf /app/data/index.db

bench:
	uv run python -m benchmarks.benchmark_basic
