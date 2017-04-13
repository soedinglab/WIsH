FROM ubuntu:14.04

RUN apt-get update
RUN apt-get install -y\
	g++\
	cmake\
	make

RUN mkdir -p /WiSH/ /data
COPY CMakeLists.txt *.h *.cpp *.in *.md *.tsv benchmark /WiSH/
RUN cd /WiSH/ && cmake . && make

WORKDIR /data
ENTRYPOINT ["/WiSH/WIsH"]
