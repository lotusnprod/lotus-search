import pytest

from api.api import autocomplete_taxa, depiction_structure
from api.models import AutocompleteTaxa, DepictionStructure
from tests.common import data_model


@pytest.mark.usefixtures("data_model")
class TestApiDepictionStructure:
    @pytest.mark.asyncio
    async def test_depiction_structure(self, data_model):
        result = await depiction_structure(DepictionStructure(structure="C"))
        assert "svg version='1.1' baseProfile='full'" in result["svg"]
        assert "FF0000" not in result["svg"]

    @pytest.mark.asyncio
    async def test_depiction_structure_highlight(self, data_model):
        result = await depiction_structure(
            DepictionStructure(structure="CO", highlight="C")
        )
        assert "svg version='1.1' baseProfile='full'" in result["svg"]
        assert "FF0000" in result["svg"]
