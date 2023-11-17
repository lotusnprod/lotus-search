from rest_framework import serializers

from .models import Reference


class ReferenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reference
        fields = ["wid", "created", "data"]
