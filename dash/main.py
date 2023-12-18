from fastapi.middleware.wsgi import WSGIMiddleware

from api.api import app
from dash_app import app as dashboard1

app.mount("", WSGIMiddleware(dashboard1.server))
