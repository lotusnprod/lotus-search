import pytest

from api.api import autocomplete_taxa
from api.models import AutocompleteTaxa
from tests.common import data_model


@pytest.mark.usefixtures("data_model")
class TestApiAutocomplete:
    @pytest.mark.asyncio
    async def test_taxa_autocomplete(self, data_model):
        result = await autocomplete_taxa(AutocompleteTaxa(taxon_name="Taxon p"), data_model)
        assert result == {"Taxon parent": 5}
