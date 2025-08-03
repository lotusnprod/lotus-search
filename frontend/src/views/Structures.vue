<template>
  <v-app>
    <v-main>
      Structures

      <h1>{{ msg }}</h1>

      <v-text-field v-model="structure"></v-text-field>
      <v-checkbox v-model="substructure_search" label="Substructure search"></v-checkbox>
      <h1>Made By Getters</h1>
      <div v-for="(structure, key) in getStructures" :key="key">
        {{key }} - {{ structure }}
      </div>
      <h1>Made By Actions</h1>
      <div v-for="(structure, key) in structures" :key="key">
        {{key}} - {{ structure }}
      </div>
    </v-main>
  </v-app>
</template>
<script setup lang="ts">
import { ref, onMounted, computed, watch } from "vue";
import { useStructureStore } from "../store/structures.ts";
const store = useStructureStore();
const structure = ref("");
const substructure_search = ref(true);
const msg = ref("Nice structures");
const getStructures = computed(() => {
  return store.getStructures;
});
const structures = computed(() => {
  return store.structures;
});
onMounted(() => {
  if (structure.value!="")
    store.fetchStructures(structure.value, substructure_search.value);
});
watch(structure, (newStructure) => {
  if (newStructure){
    store.fetchStructures(newStructure, substructure_search.value);
  }
});
</script>
