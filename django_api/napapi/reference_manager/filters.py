from django_filters import rest_framework as filters

from .models import Reference


class ReferenceFilter(filters.FilterSet):
    wid = filters.Filter(lookup_expr='exact')

    class Meta:
        model = Reference
        fields = ['wid']
