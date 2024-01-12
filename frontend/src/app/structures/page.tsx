'use client'

import StructureSearch from "@/components/StructureSearch";
import StructureResult from "@/components/StructuresResult";
import {useState} from "react";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";

export default function Page() {
    const [searchQuery, setSearchQuery] = useState<StructureSearchQuery>({smiles: 'c1ccccc1', substructureSearch: false});
    return (
        <main className="flex min-h-screen flex-col items-center justify-between p-24">
            <div className="z-10 max-w-5xl w-full items-center justify-between font-mono text-sm lg:flex">
                <StructureSearch onSearchSubmit={setSearchQuery} />
            </div>
            <div className="z-10 max-w-5xl w-full items-center justify-between font-mono text-sm lg:flex">
                <StructureResult searchQuery={searchQuery} />
            </div>
        </main>
    )
}
