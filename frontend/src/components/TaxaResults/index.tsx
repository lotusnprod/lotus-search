'use client'

import React, {useEffect, useState} from 'react';
import {fetchStructures, fetchTaxa} from "@/services/apiService";
import {components} from "@/interfaces/lotus_api";
import Structure from "@/components/Structure";
import Pagination from "@/components/Pagination";
import {Grid, Sheet} from "@mui/joy";
import {LotusAPIItem} from "@/interfaces/schemas";

interface TaxonResultQuery {
    searchQuery: LotusAPIItem;
}


const ITEMS_PER_PAGE = 20; // Set the number of items per page


const TaxaResults: React.FC<TaxonResultQuery> = ({searchQuery}) => {
    const [currentPage, setCurrentPage] = useState(1);
    const [apiData, setApiData] = useState<components["schemas"]["TaxonResult"] | null>(null);
    const [loading, setLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    useEffect(() => {
        if (searchQuery.taxon_name == "") return
        setLoading(true)
        fetchTaxa(searchQuery).then(setApiData)
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
    const totalPages = apiData && apiData.taxa ? Math.ceil(Object.keys(apiData.taxa).length / ITEMS_PER_PAGE) : 0;

    // Get the current items
    const currentItems: Taxa[] = apiData && apiData.taxa
        ? Object.entries(apiData.taxa)
            .slice((currentPage - 1) * ITEMS_PER_PAGE, currentPage * ITEMS_PER_PAGE)
        : [];
    if (currentItems.length == 0) {
        return <div>No results</div>
    }
    return (<Sheet>
        Taxon searching
        <Grid container spacing={2} sx={{flexGrow: 1}}>
            {currentItems.map(([index, item]) => (
                <div key={index}>{ item.name }</div>
            ))}
        </Grid>
        <Pagination currentPage={currentPage} totalPages={totalPages} onPageChange={setCurrentPage}/>
    </Sheet>)
}

export default TaxaResults
