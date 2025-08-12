from dash_app import app as dashboard1
from fastapi.middleware.wsgi import WSGIMiddleware

from api.api import app

app.mount("", WSGIMiddleware(dashboard1.server))
