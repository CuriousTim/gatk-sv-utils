docker_name := gatk-sv-utils
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

docker_tag := $(if $(DOCKER_REPO),$(DOCKER_REPO)/)$(docker_name):$(date)-$(branch)-$(commit_sha)

.PHONY: all
all: ;

.PHONY: docker-build
docker-build:
	docker build -t $(docker_tag) .

.PHONY: docker-push
docker-push:
	docker push $(docker_tag)
