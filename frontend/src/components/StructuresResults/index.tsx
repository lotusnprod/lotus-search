'use client'

import React, {useEffect, useState} from 'react';
import {fetchStructures} from "@/services/apiService";
import {components} from "@/interfaces/lotus_api";
import Structure from "@/components/Structure";
import Pagination from "@/components/Pagination";
import {Grid, Sheet} from "@mui/joy";
import {LotusAPIItem} from "@/interfaces/schemas";

interface StructureResultProps {
    searchQuery: LotusAPIItem;
}


const ITEMS_PER_PAGE = 20; // Set the number of items per page


const StructureResult: React.FC<StructureResultProps> = ({searchQuery}) => {
    const [currentPage, setCurrentPage] = useState(1);
    const [apiData, setApiData] = useState<components["schemas"]["StructureResult"] | null>(null);
    const [loading, setLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    useEffect(() => {
        if (searchQuery.structure == "") return
        setLoading(true)
        fetchStructures(searchQuery).then(setApiData)
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
    if (currentItems.length == 0) {
        return <div>No results</div>
    }
    return (<Sheet>
        Structure searching
        <Grid container spacing={2} sx={{flexGrow: 1}}>
            {currentItems.map(([index, structure]) => (
                structure?.smiles ? <Grid key={"grid_structure_" + index} xs={4}>
                        <Structure key={"structure_" + index} id={index} structure={structure.smiles}
                                   highlight={searchQuery.smiles}/>
                    </Grid>
                    : null
            ))}
        </Grid>
        <Pagination currentPage={currentPage} totalPages={totalPages} onPageChange={setCurrentPage}/>
    </Sheet>)
}

export default StructureResult
