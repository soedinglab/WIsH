FROM ubuntu:14.04

RUN apt-get update
RUN apt-get install -y\
	g++\
	cmake\
	make

RUN mkdir -p /WiSH/build /data
COPY CMakeLists.txt *.h *.cpp *.in *.md *.tsv benchmark /WiSH/
RUN cd /WiSH/build && cmake .. && make

WORKDIR /data
ENTRYPOINT ["/WiSH/build/WIsH"]
