FROM python:3.13.5-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxi6 \
    libxrender1 \
    libxtst6 \
    && rm -rf /var/lib/apt/lists/*
RUN pip install uv
COPY uv.lock pyproject.toml ./
RUN --mount=type=cache,target=/root/.cache/pdm uv sync
RUN mkdir /app
RUN adduser --system --no-create-home nonroot
USER nonroot
WORKDIR /app
CMD ["gunicorn", "--workers=4", "--threads=1", "-b", "0.0.0.0:8502", "main:server"]
