FROM public.ecr.aws/docker/library/node:20
RUN npm i -g @redocly/cli@1.6.0 openapi-typescript@6.7.3
RUN mkdir -p /spec
WORKDIR /spec
RUN chown -R node:node /spec
USER node
CMD ["bash"]
