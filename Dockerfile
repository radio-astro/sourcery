FROM radioastro/base:0.2

MAINTAINER gijsmolenaar@gmail.com

RUN apt-get update && \
    apt-get install -y \
        lofar \
        python-numpy \
        python-scipy \
        python-astlib \
        python-tigger \
    &&  \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ADD . /tmp/sourcery

RUN cd /tmp/sourcery && pip install .

CMD python /usr/local/bin/sourcery
