from fastapi.middleware.wsgi import WSGIMiddleware

from api import app as api
from dash_app import app as dashboard1

app = api

app.mount("", WSGIMiddleware(dashboard1.server))
