// Composables
import {createRouter, createWebHistory} from 'vue-router'

const routes = [
    {
        path: '/',
        component: () => import('@/layouts/default/Default.vue'),
        children: [
            {
                path: '',
                name: 'Home',
                component: () => import(/* webpackChunkName: "home" */ '@/views/Home.vue'),
            },
            {
                path: '/structures',
                name: 'Structures',
                component: () => import(/* webpackChunkName: "structures" */ '@/views/Structures.vue'),
            },
            {
                path: '/taxa',
                name: 'Taxa',
                component: () => import(/* webpackChunkName: "taxa" */ '@/views/Taxa.vue'),
            }
        ],
    },
]

const router = createRouter({
    history: createWebHistory(process.env.BASE_URL),
    routes,
})

export default router
