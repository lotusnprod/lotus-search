'use client'

import React, {useRef, useState} from "react";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";
import KetcherLotus, {KetcherLotusMethods} from "@/components/KetcherLotus";
import dynamic from "next/dynamic";

interface StructureSearchProps {
    onSearchSubmit: (searchValue: StructureSearchQuery) => void;
}

const StructureSearch: React.FC<StructureSearchProps> =
    ({onSearchSubmit}) => {
        const [substructureSearch, setSubstructureSearch] = useState(false);
        const ketcherRef = useRef<KetcherLotusMethods>(null);
        const handleSubmit = async (event: React.FormEvent<HTMLFormElement>) => {
            event.preventDefault();
            if (ketcherRef.current) {
                const smiles = await ketcherRef.current.getValue()
                if (smiles != "")
                    onSearchSubmit({smiles: smiles, substructureSearch: substructureSearch})
            } else {
                console.log("Nothing to search for")
            }
        };

        return (
            <main>
                <div className="w-full">
                    <KetcherLotus ref={ketcherRef}/>
                </div>
                <form className="flex flex-col space-y-4" onSubmit={handleSubmit}>
                    <label className="inline-flex items-center">
                        <input className='bg-gray-200 text-black shadow-inner rounded-l p-2 flex-1' id='smiles'
                               type='checkbox'
                               aria-label='Substructure search'
                               onChange={e => setSubstructureSearch(!substructureSearch)}
                        />
                        <span className="ml-2 text-gray-700">Substructure search</span>
                    </label>

                    <button className='bg-blue-600 hover:bg-blue-700 duration-300 text-white shadow p-2 rounded-r'
                            type='submit'>
                        Search structure
                    </button>
                    {}
                </form>
            </main>
        )
    }
export default dynamic(() => Promise.resolve(StructureSearch), {
    ssr: false
})
