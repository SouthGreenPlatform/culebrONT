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
	@echo "CulebrONT installed in ${TARGET_DIR}/bin/canu"

.PHONY=build
build:
ifneq ($(shell id -u), 0)
	@echo "You must be root to perform this action.";
else
	@echo "Build all singularity"
	cd Containers; \
	for def in `ls *.def`; do \
		simg=$${def%%.def}.sif;\
		if test -f $${simg}; then \
			echo "Skip build $${simg} already exist !!"; \
		else \
			echo "Build $${simg} from def: $${def}";\
			singularity build $${simg} $${def};\
		fi; \
	done
	@echo "Build all conda"
	for yaml in `ls ./envs/*.yaml`; do \
		env_name=$${yaml%%.yaml};\
		if test -d ./envs/$${env_name}; then \
			echo "Skip build ./envs/$${env_name} already exist !!"; \
		else \
			echo "Build ./envs/$${env_name} from yaml: $${yaml}";\
			conda env create -f $${yaml} -p ./envs/$${env_name};\
		fi; \
	done
endif

.PHONY=install
install:
	@echo "Try to install";

	@echo "Build all conda";
	cd envs; \
	for yaml in `ls *.yaml`; do \
		env_name=$${yaml%%.yaml};\
		if test -d $${env_name}; then \
			echo "Skip build $${env_name} already exist !!"; \
		else \
			echo "Build $${env_name} from yaml: $${yaml}";\
			conda env create -f $${yaml} -p $${env_name};\
			echo "Build DONE";\
		fi; \
	done
