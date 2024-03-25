import pytest
from rdkit.Chem import Descriptors

from api.api import get_descriptors


class TestApiDescriptors:
    @pytest.mark.asyncio
    async def test_get_descriptors(self):
        result = await get_descriptors()
        assert result == [desc[0] for desc in Descriptors._descList]
