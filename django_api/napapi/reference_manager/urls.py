from django.urls import path

from .views import (ReferenceCreate, ReferenceDelete, ReferenceDetail,
                    ReferenceDetailWid, ReferenceList, ReferenceUpdate)

urlpatterns = [
    path("create/", ReferenceCreate.as_view(), name="create-reference"),
    path("", ReferenceList.as_view()),
    path("<int:pk>/", ReferenceDetail.as_view(), name="retrieve-reference"),
    path(
        "wid/<int:wid>/", ReferenceDetailWid.as_view(), name="retrieve-reference-by-wid"
    ),
    path("update/<int:pk>/", ReferenceUpdate.as_view(), name="update-reference"),
    path("delete/<int:pk>/", ReferenceDelete.as_view(), name="delete-reference"),
]
