'use client'

import React, {useEffect, useState} from 'react';
import {fetchStructures} from "@/services/apiService";
import {components} from "@/interfaces/lotus_api";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";
import Structure from "@/components/Structure";
import Pagination from "@/components/Pagination";

interface StructureResultProps {
    searchQuery: StructureSearchQuery;
}


const ITEMS_PER_PAGE = 20; // Set the number of items per page


const StructureResult: React.FC<StructureResultProps> = ({searchQuery}) => {
    const [currentPage, setCurrentPage] = useState(1);
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

     // Calculate the total number of pages
  const totalPages = apiData && apiData.structures ? Math.ceil(Object.keys(apiData.structures).length / ITEMS_PER_PAGE) : 0;

  // Get the current items
  const currentItems = apiData && apiData.structures
    ? Object.entries(apiData.structures)
        .slice((currentPage - 1) * ITEMS_PER_PAGE, currentPage * ITEMS_PER_PAGE)
    : [];

    return (<div>Structure searching
     <div className="flex flex-wrap">
        {currentItems.map(([index, structure]) => (
          structure?.smiles ? <Structure key={"structure_" + index} id={index} structure={structure.smiles} highlight={searchQuery.smiles}/> : null
        ))}
      </div>
      <Pagination currentPage={currentPage} totalPages={totalPages} onPageChange={setCurrentPage} />

</div>)
}

export default StructureResult
