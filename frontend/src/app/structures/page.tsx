'use client'

import StructureSearch from "@/components/StructureSearch";
import StructureResults from "@/components/StructuresResults";
import {useState} from "react";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";
import {LotusAPIItem} from "@/interfaces/schemas";

export default function Page() {
    const [searchQuery, setSearchQuery] = useState<LotusAPIItem>({structure: '', substructure_search: false});
    return (
        <main className="flex min-h-screen flex-col items-center justify-between p-24">
            <div className="z-10 max-w-5xl w-full items-center justify-between font-mono text-sm">
                <StructureSearch onSearchSubmit={setSearchQuery} />
            </div>
            <div className="z-10 max-w-5xl w-full items-center justify-between font-mono text-sm">
                <StructureResults searchQuery={searchQuery} />
            </div>
        </main>
    )
}