import {components} from "@/interfaces/lotus_api";

export const fetchStructures = async (input: components["schemas"]["item"]): Promise<components["schemas"]["StructureResult"]> => {
   const response = await fetch('/api/v1_0/structures', {
        method: 'POST',
        headers: {
             'Content-Type': 'application/json'
        },
        body: JSON.stringify(input)
     });

    if (!response.ok) {
        throw new Error(response.statusText);
    }
    return await response.json();
}
