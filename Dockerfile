FROM docker.io/library/python:3.11-slim
RUN apt update && apt install -y libxrender1 libxtst6 libxi6
RUN pip install poetry
COPY poetry.lock pyproject.toml ./
RUN poetry config virtualenvs.create false \
    && poetry install --no-interaction --no-ansi
RUN mkdir /app
RUN adduser --system --no-create-home nonroot
USER nonroot
WORKDIR /app
CMD [ "gunicorn", "--workers=2", "--threads=1", "-b 0.0.0.0:8502", "main:server"]
