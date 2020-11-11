lat = [48.8708 51.5188 41.9260 40.4312 52.523 37.982]; 
lon = [2.4131 -0.1300 12.4951 -3.6788 13.415 23.715]; 


plot(lon,lat,'.r','MarkerSize',20) 
plot_google_map( 'maptype', 'satellite')