import json

from django.db import models
import requests
from django.conf import settings



class Reference(models.Model):
    wid = models.IntegerField(unique=True, null=True, blank=True, db_index=True)
    created = models.DateField(auto_now_add=True)
    data = models.JSONField()

    def title(self):
        if "title" in self.data and self.data["title"] is not None:
            return self.data["title"]
        return "Untitled"

    @staticmethod
    def gather_from_wikidata(wid: int):
        query = """
            SELECT ?p ?o WHERE {
            VALUES ?article_type { wd:Q580922 wd:Q13442814 }
              wd:Q%WID% wdt:P31 ?article_type ; 
                         ?p ?o.
            }
            """.replace("%WID%", str(wid))
        result = json.loads(requests.get(settings.WD_URL,
                            params={'query': query},
                            headers={'Accept': 'application/json'}).text)

        data = {}

        if len(result["results"]["bindings"]) < 2:
            raise Reference.DoesNotExist("Impossible to find the reference in Wikidata")

        for binding in result["results"]["bindings"]:
            if binding["p"]["value"] not in data:
                data[binding["p"]["value"]] = []
            data[binding["p"]["value"]].append(binding["o"]["value"])

        ref = Reference(wid=wid, data=data)
        ref.save()
        return ref

    def __str__(self):
        if self.data is None:
            if self.wid is not None:
                return f"Q{self.wid}"

        return self.title()
