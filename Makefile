

.PHONY: start stop tests fix

start:
	docker compose up --build -d

stop:
	docker compose down


tests: start
	docker compose run backend poetry run pytest --cov-config=.coveragerc --cov=. -n auto

fix: start
	docker compose run backend poetry run black .
	docker compose run backend poetry run isort .

.PHONY: clean_data clean_db
clean_data:
	docker compose run backend rm -rf /app/data/*

clean_db:
	docker compose run backend rm -rf /app/data/index.db
