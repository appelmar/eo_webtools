FROM ubuntu:xenial
MAINTAINER Marius Appel <marius.appel@uni-muenster.de>

RUN apt-get update && apt-get install -y software-properties-common libboost-all-dev cmake g++ golang-go git supervisor nano
RUN add-apt-repository -y ppa:ubuntugis/ppa && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 314DF160 && apt-get update && apt-get install -y libgdal-dev gdal-bin libproj-dev python-gdal


WORKDIR /opt
RUN git clone https://github.com/appelmar/eo_webtools

WORKDIR /opt/eo_webtools/eo_webtile
RUN mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ../ && make -j 2 && chmod +x eo_webtile && cp -l eo_webtile /usr/bin/

WORKDIR /opt/eo_webtools/s2vrt
RUN mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ../ && make -j 2 && chmod +x s2vrt && cp -l s2vrt /usr/bin/

ENV GOPATH=/opt/go
WORKDIR /opt/eo_webtools/sentinel_webview
RUN go get "github.com/gorilla/mux"
RUN go build 

COPY s2scan_dir.py /opt/s2scan_dir.py  
RUN mkdir /opt/VRT 

COPY supervisord.conf /opt/
EXPOSE 8080

CMD ["/usr/bin/supervisord", "-c", "/opt/supervisord.conf"]


# docker build -t "appelmar/sentinel_webview" .
# run with docker run -p 8080:8080 -v ~/eodata/Sentinel2/:/data   appelmar/sentinel_webview


