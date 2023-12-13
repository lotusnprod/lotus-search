FROM node:21.4
ENV PNPM_HOME="/pnpm"
ENV PATH="$PNPM_HOME:$PATH"
RUN corepack enable
WORKDIR /app
COPY package.json .
RUN pnpm install
COPY tsconfig.json .
COPY tsconfig.node.json .
COPY vite.config.ts .
COPY package.json .
COPY pnpm-lock.yaml .
COPY index.html .
CMD ["pnpm", "run", "dev"]
