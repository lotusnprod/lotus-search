.PHONY: all clean test api

# Standard required targets
all: api

clean:
	@echo "Cleaning API documentation..."
	@$(MAKE) -C api clean 2>/dev/null || echo "No clean target in api subdirectory"

test: api
	@echo "API documentation generated successfully"

# Your existing target
api:
	@echo "Generating API documentation..."
	@$(MAKE) -C api
