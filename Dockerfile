FROM  mengchen18/xcmsviewer_shiny:0.6.9

RUN mkdir /home/shiny/rtkiapp 

COPY ./ /home/shiny/rtkiapp/

RUN R -e 'install.packages(c("heatmaply", "randomcoloR"))'

CMD ["R", "-e", "shiny::runApp('/home/shiny/rtkiapp', host = '0.0.0.0', port = 3838)"]
