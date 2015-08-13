library(shiny)
library(leaflet)
# Calculates the geodesic distance between two points specified by radian latitude/longitude using
# Vincenty inverse formula for ellipsoids (vif)
# Credit to Mario Pineda-Krch and Chris Veness
# http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.vif <- function(long1, lat1, long2, lat2) {
        
        # WGS-84 ellipsoid parameters
        a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
        b <- 6356752.314245  # ength of minor axis of the ellipsoid (radius at the poles)
        f <- 1/298.257223563 # flattening of the ellipsoid
        
        L <- long2-long1 # difference in longitude
        U1 <- atan((1-f) * tan(lat1)) # reduced latitude
        U2 <- atan((1-f) * tan(lat2)) # reduced latitude
        sinU1 <- sin(U1)
        cosU1 <- cos(U1)
        sinU2 <- sin(U2)
        cosU2 <- cos(U2)
        
        cosSqAlpha <- NULL
        sinSigma <- NULL
        cosSigma <- NULL
        cos2SigmaM <- NULL
        sigma <- NULL
        
        lambda <- L
        lambdaP <- 0
        iterLimit <- 100
        while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
                sinLambda <- sin(lambda)
                cosLambda <- cos(lambda)
                sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                                          (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
                if (sinSigma==0) return(0)  # Co-incident points
                cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
                sigma <- atan2(sinSigma, cosSigma)
                sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
                cosSqAlpha <- 1 - sinAlpha*sinAlpha
                cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
                if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
                C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
                lambdaP <- lambda
                lambda <- L + (1-C) * f * sinAlpha *
                        (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
                iterLimit <- iterLimit - 1
        }
        if (iterLimit==0) return(NA)  # formula failed to converge
        uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
        A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
        B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
        deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                                         B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
        s <- b*A*(sigma-deltaSigma) / 1000
        
        return(s) # Distance in km
}

#plotmap = function(lon1, lat1, lon2, lat2){
#         map3 <- Leaflet$new()
#         # Convert Radian to Degree
#         lon1=lon1/pi *180
#         lon2=lon2/pi *180
#         lat1=lat1/pi *180
#         lat2=lat2/pi *180
#         
#         map3$setView(c(lat1, lon1), zoom= 2) 
#         map3$marker(c(lat1, lon1), bindPopup = "<p> Origin </p>") 
#         map3$marker(c(lat2, lon2), bindPopup = "<p> Destination </p>") 
#         map3
        #map3$save('fig/map3.html ', cdn = TRUE) 
        #cat( '<iframe src="fig/map3.html" width=100%, height=600></iframe>') 
        
#}

shinyServer(
        function(input, output) {
                output$gcMap <- renderLeaflet({
                        #x=c(input$lat1, input$lat2)
                        #y=c(input$lon1, input$lon2)
                        #plot(x,y)
                        #lines(x,y)
                        dist=gcd.vif(input$lon1, input$lat1, input$lon2, input$lat2)
                        # Convert Radian to Degree
                        lon1=input$lon1/pi *180
                        lon2=input$lon2/pi *180
                        lat1=input$lat1/pi *180
                        lat2=input$lat2/pi *180
                        # Render leaflet map
                        map3 <- leaflet() %>% addTiles() %>% addMarkers(lng=lon1, lat=lat1, popup="Origin") %>% addMarkers(lng=lon2, lat=lat2, popup="Destination") %>% setView(lng=mean(c(lon1, lon2)), lat=mean(c(lat1, lat2)), zoom = 2)
                })
                
                output$title= renderPrint(paste("GREAT CIRCLE DIST (km)", gcd.vif(input$lon1, input$lat1, input$lon2, input$lat2)))
        }
)