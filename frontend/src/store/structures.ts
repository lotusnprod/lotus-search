import { defineStore } from 'pinia'
import axios from "axios"

export const useStructureStore = defineStore("structure", {
    state: () => ({
        structures: [],
    }),
    getters: {
        getStructures(state){
            return state.structures
        }
    },
    actions: {
        async fetchStructures(structure: string, substructure_search: boolean=true) {
            try {
                const data = await axios.post('/api/v1_0/structures',
                    {"molecule": structure, "substructure_search": substructure_search})
                this.structures = data.data.structures
            }
            catch (error) {
                console.log(error)
            }
        }
    },
})
