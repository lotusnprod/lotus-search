const CopyPlugin = require("copy-webpack-plugin");
/** @type {import('next').NextConfig} */
const nextConfig = {
    webpack(config, {isServer}) {
        config.plugins.push(
            new CopyPlugin({
                patterns: [
                    {
                        from: "node_modules/@rdkit/rdkit/dist/RDKit_minimal.wasm",
                        to: "static/chunks"
                    }
                ]
            })
        );

        if (!isServer) {
            config.resolve.fallback = {
                fs: false
            };
        }

        return config;
    },


    async rewrites() {
        return [
            {
                source: '/api/:path*',
                destination: 'http://backend:5000/:path*' // Proxy to Backend
            }
        ]
    }
}

module.exports = nextConfig
