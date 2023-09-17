from django.contrib import admin

from .models import Reference

class ReferenceAdmin(admin.ModelAdmin):
    list_display = ["wid", "title"]

admin.site.register(Reference, ReferenceAdmin)
