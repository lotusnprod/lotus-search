<template>
  <v-app-bar title="LOTUS"></v-app-bar>

  <v-navigation-drawer>
    <v-list>
      <v-list-item :to='{name: "Home"}'>Home</v-list-item>
      <v-list-item to="/structures">Search structures</v-list-item>
      <v-list-item to="/taxa">Search taxa</v-list-item>
    </v-list>
  </v-navigation-drawer>

  <v-main class="d-flex align-center justify-center" style="min-height: 300px;">
    <default-view />
  </v-main>
</template>

<script lang="ts" setup>
  import DefaultView from './View.vue'
</script>
