

.PHONY: start tests

start:
	docker compose up --build -d

tests: start
	docker compose run backend poetry run pytest --cov-config=.coveragerc --cov=. -n auto

fix: start
	docker compose run backend poetry run black .
	docker compose run backend poetry run isort .
