<script setup lang="ts">
import { ref } from 'vue'

defineProps<{ msg: string }>()

const count = ref(0)
</script>

<template>
  <h1>{{ msg }}</h1>

  <v-btn @click="count++">count is {{ count }}</v-btn>
</template>

<style scoped>

</style>
