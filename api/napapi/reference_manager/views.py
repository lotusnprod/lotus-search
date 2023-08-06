from django.shortcuts import get_object_or_404, render
from rest_framework import generics
from django_filters.rest_framework import DjangoFilterBackend

from .filters import ReferenceFilter
from .models import Reference
from .serializers import ReferenceSerializer


class ReferenceCreate(generics.CreateAPIView):
    queryset = Reference.objects.all(),
    serializer_class = ReferenceSerializer


class ReferenceList(generics.ListAPIView):
    queryset = Reference.objects.all()
    serializer_class = ReferenceSerializer
    filter_backends = (DjangoFilterBackend,)
    filterset_class = ReferenceFilter

class ReferenceDetail(generics.RetrieveAPIView):
    queryset = Reference.objects.all()
    serializer_class = ReferenceSerializer

class ReferenceDetailWid(generics.RetrieveAPIView):
    queryset = Reference.objects.all()
    serializer_class = ReferenceSerializer

    def get_object(self):
        wid = self.kwargs['wid']
        try:
            obj = Reference.objects.get(wid=wid)
        except Reference.DoesNotExist:
            obj = Reference.gather_from_wikidata(wid)

        return obj

class ReferenceUpdate(generics.RetrieveUpdateAPIView):
    queryset = Reference.objects.all()
    serializer_class = ReferenceSerializer

class ReferenceDelete(generics.RetrieveDestroyAPIView):
    queryset = Reference.objects.all()
    serializer_class = ReferenceSerializer

