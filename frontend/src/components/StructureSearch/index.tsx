'use client'

import {useState} from "react";
import {StructureSearchQuery} from "@/interfaces/structure_search_query";

interface StructureSearchProps {
    onSearchSubmit: (searchValue: StructureSearchQuery) => void;
}

const StructureSearch: React.FC<StructureSearchProps> =
    ({onSearchSubmit}) => {
        const [inputValue, setInputValue] = useState('c1ccccc1');
        const [substructureSearch, setSubstructureSearch] = useState(false);

        const handleSubmit = (event: React.FormEvent<HTMLFormElement>) => {
            event.preventDefault();
            onSearchSubmit({smiles: inputValue, substructureSearch: substructureSearch});
        };

        return (
            <main>
                <form className="flex flex-col space-y-4" onSubmit={handleSubmit}>
                    <input className='bg-gray-200 text-black shadow-inner rounded-l p-2 flex-1' id='smiles' type='text'
                           aria-label='structure smiles' placeholder='SMILES string.'
                           value={inputValue}
                           onChange={(e) => setInputValue(e.target.value)}
                    />
                    <label className="inline-flex items-center">
                        <input className='bg-gray-200 text-black shadow-inner rounded-l p-2 flex-1' id='smiles'
                               type='checkbox'
                               aria-label='Substructure search'
                               value={inputValue}
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

export default StructureSearch;