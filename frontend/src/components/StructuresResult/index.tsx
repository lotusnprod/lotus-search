'use client'

import React, {useEffect, useState} from 'react';
import {fetchStructures} from "@/services/apiService";
import {components} from "@/interfaces/lotus_api";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";
import Structure from "@/components/Structure";

interface StructureResultProps {
    searchQuery: StructureSearchQuery;
}

const StructureResult: React.FC<StructureResultProps> = ({searchQuery}) => {
    const [apiData, setApiData] = useState<components["schemas"]["StructureResult"] | null>(null);
    const [loading, setLoading] = useState<boolean>(true);
    const [error, setError] = useState<string | null>(null);
    useEffect(() => {
        if (searchQuery.smiles == "") return
        fetchStructures({
            "structure": searchQuery.smiles,
            "substructure_search": searchQuery.substructureSearch
        }).then(setApiData)
            .catch((error) => setError(error.message))
            .finally(() => {
                    setLoading(false)
                    setError(null)
                }
            );
    }, [searchQuery]);

    if (loading) return <div>Loading...</div>;
    if (error) return <div>Error: {error}</div>;

    return (<div>Structure searching
    <div  className="flex flex-wrap">
        {apiData && apiData.structures && Object.entries(apiData.structures).map(([index, structure]) => (
            structure?.smiles ? <Structure id={index} structure={structure.smiles} highlight={searchQuery.smiles}/> : ''
))}
        </div>

{
    apiData && <pre>{JSON.stringify(apiData, null, 2)}</pre>
}
</div>)
}

export default StructureResult
