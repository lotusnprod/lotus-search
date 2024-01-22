'use client'

import React, {useState} from "react";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";
import AutocompleteTaxa from "@/components/AutocompleteTaxa";

export default function Page() {
    const [searchQuery, setSearchQuery] = useState<StructureSearchQuery>({smiles: '', substructureSearch: false});
    return (
        <main className="flex min-h-screen flex-col items-center justify-between p-24">
            <div className="z-10 max-w-5xl w-full items-center justify-between font-mono text-sm">
                <AutocompleteTaxa />
            </div>
            <div className="z-10 max-w-5xl w-full items-center justify-between font-mono text-sm">
                The results
            </div>
        </main>
    )
}
