DOCKER_REPO=radioastro/sourcery:1.2.5

.PHONY: build clean

all: build

build:
	docker build -t ${DOCKER_REPO} .

run:
	docker run -ti ${DOCKER_REPO}

clean:
	docker rmi ${DOCKER_REPO}

upload: build
	docker push ${DOCKER_REPO}
