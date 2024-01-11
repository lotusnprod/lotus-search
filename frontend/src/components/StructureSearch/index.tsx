export default function StructureSearch() {
    return(
    <main>
        <form className='flex'>
            <input className='bg-gray-200 text-black shadow-inner rounded-l p-2 flex-1' id='smiles' type='text'
                   aria-label='structure smiles' placeholder='SMILES string.'/>
            <button className='bg-blue-600 hover:bg-blue-700 duration-300 text-white shadow p-2 rounded-r'
                    type='submit'>
                Search structure
            </button>
        </form>
    </main>
    )
}
