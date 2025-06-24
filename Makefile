project_name := gatk-sv-utils
commit_sha := $(shell git rev-parse --short HEAD)
ifneq ($(.SHELLSTATUS), 0)
$(error failed to get git commit hash)
endif

date := $(shell date --utc +%Y%m%d)
ifneq ($(.SHELLSTATUS), 0)
$(error failed to get date)
endif

branch := $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(.SHELLSTATUS), 0)
$(error failed to get git branch)
endif

docker_tag = $(if $(DOCKER_REPO),$(DOCKER_REPO)/)$(project_name):$(date)-$(docker_name)-$(commit_sha)

define build-docker =
docker build -t $(docker_tag) -f $< .
endef

define push-docker =
docker push $(docker_tag)
endef

.PHONY: all
all: ;

.PHONY: docker
docker-build: docker-base docker-r

.PHONY: docker-base
docker-base: docker_name := base
docker-base: docker/base/Dockerfile
	$(build-docker)
ifneq ($(strip $(DOCKER_REPO)),)
	$(push-docker)
endif

.PHONY: docker-r
docker-r: docker_name := r
docker-r: docker/r/Dockerfile
	$(build-docker)
ifneq ($(strip $(DOCKER_REPO)),)
	$(push-docker)
endif
