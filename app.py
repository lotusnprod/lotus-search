"""Application entrypoint.

Provides a `main()` function so the pyproject script `start = "app:main"` works.
It simply runs the FastAPI application defined in `api.api` using uvicorn.
(No change to existing API or Dash functionality.)
"""
from __future__ import annotations

import uvicorn


def main() -> None:  # pragma: no cover - thin wrapper
    """Run the FastAPI (and mounted Dash) application.

    This is intentionally minimal to avoid side effects and keep
    backwards compatibility. All application construction happens
    in `api.api`.
    """
    uvicorn.run(
        "api.api:app",
        host="0.0.0.0",
        port=8000,
        reload=False,
        factory=False,
    )


if __name__ == "__main__":  # pragma: no cover
    main()

