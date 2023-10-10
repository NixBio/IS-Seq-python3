# docker build --network 'host' --no-cache -f $HOME/Aimin/IS-Seq-python3/Docker/Dockerfile -t aiminy/isseq:1.8 .
# to RUN
# docker run -it aiminy/isseq:2.1 bash

# docker build --no-cache -f $HOME/Aimin/IS-Seq-python3/Docker/Dockerfile -t aiminy/isseq:2.1 .

# docker run --mount type=bind,source=/local_scratch,target=/home -it aiminy/isseq:2.1 bash

# docker system prune --all --force --volumes

From aiminy/isseq:2.1

MAINTAINER Aimin Yan <aimin.at.work@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

COPY /usr/src/IS-Seq-python3/R

