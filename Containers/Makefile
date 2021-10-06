help:
	@echo "available subcommands:"
	@echo "\t build: use to build all sif from .def file"
	@echo "\t install: use to autoconfigure .sif to CulebrONT  cluster-config profile"
	@echo "\t all: run build an install"

.PHONY=all
all:
	make build
	make install
	@echo ""
	@echo "Success!"
	@echo "CulebrONT installed in ${TARGET_DIR}"

.PHONY=build
build:
ifneq ($(shell id -u), 0)
	@echo "You must be root to perform this action.";
else
	@echo "Build all singularity"
	for def in `ls *.def`; do \
		simg=$${def%%.def}.sif;\
		if test -f $${simg}; then \
			echo "Skip build $${simg} already exist !!"; \
		else \
			echo "Build $${simg} from def: $${def}";\
			singularity build $${simg} $${def};\
		fi; \
	done
endif

.PHONY=install
install:
	@echo "Try to install";