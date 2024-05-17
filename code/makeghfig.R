# dependencies
# library(LambertW)
# library(gsl)
# library(nleqslv)

# Code from Yuan Yan <yuan.yan@dal.ca>
# Reference:
# Yan, Yuan, and Marc G. Genton. ‘The Tukey g-and-h Distribution’. Significance 16, no. 3 (2019): 12–13. https://doi.org/10.1111/j.1740-9713.2019.01273.x.

gh_trans <- function(z,g=0.2,h=0.1){
  if(g==0)
    x<-z*exp(h*z^2/2)
  else
    x<-1/g*(exp(g*z)-1)*exp(h*(z^2)/2)
  return(x)
}

#pdf: density function
den_tukey <- function(y,g=0.2,h=0.1,xi=0,omega=1,u=0,sig=1){
  y=(y-xi)/omega  
  z=tukey_inv(y,g,h)
  #   for(i in 1:length(y)){
  #     # z[i] <- uniroot(f,c(-10,10))$root
  #     z[i]<-nleqslv(0,fn=function(x){tukey(x,g,h)-y[i]})$x
  #   }
  dnorm(z,mean=u,sd=sig)/dtukey(z,g,h)/omega  
}

#inverse TGH transformation
tukey_inv <- function(y,g=0.2,h=0.1){
  if(g==0 & h==0)
    z=y
  else if(g==0)
    z<-sign(y)*sqrt(lambert_W0(h*y^2)/h) #better & 143 times faster
  else if(h==0){
    # if(g>0)
    #   y=y[y>-1/g]
    # else
    #   y=y[y<-1/g]
    z<-log(g*y+1)/g
  }
  else{
    z=numeric(length(y))
    for(i in 1:length(y)){
      # z[i] <- uniroot(f,c(-10,10))$root
      z[i]<-nleqslv(0,fn=function(x){gh_trans(x,g,h)-y[i]})$x
    }
  }
  return(z)
}

dtukey <- function(z,g=0.2,h=0.1){
  if(g==0)
    x<-exp(h*z^2/2)*(1+h*z^2)
  else
    x<-exp(g*z+h*z^2/2)+h*z/g*(exp(g*z)-1)*exp(h*(z^2)/2)
  return(x)
}

# return pdf of g-and-h distribution in which g varies and h is constant
plot_g_pdf <- function(gvec = seq(0, 1, 0.1), x = seq(-6, 6, 0.01), h = 0){
ng <- length(gvec)
nx <- length(x) 
ghpdf <- array(NA, dim = c(nx, ng))  

for(G in 1:ng){
  ghpdf[,G] <- den_tukey(x, g=gvec[G], h=h)
}

df <- tibble(density = as.vector(ghpdf),
             x = rep(x, ng),
             g = factor(rep(gvec, each = nx)))

out <- ggplot(df, aes(x = x, y = density, colour = g)) + theme_bw() +
  geom_line(linewidth = 1) +
  scale_colour_viridis_d(end = 0.9, option = "B") +
  scale_x_continuous(breaks = seq(-4, 5, 1)) +
  coord_cartesian(xlim = c(-4, 5)) +
  labs(x = "Observations",
       y = "Density") +
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(.1, .55),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(reverse=TRUE, override.aes = list(linewidth = 3)))

out
}

# return pdf of g-and-h distribution in which h varies and g is constant
plot_h_pdf <- function(hvec = seq(0, .3, 0.05), x = seq(-6, 6, 0.01), g = 0){
  nh <- length(hvec)
  nx <- length(x) 
  ghpdf <- array(NA, dim = c(nx, nh))  
  
  for(H in 1:nh){
    ghpdf[,H] <- den_tukey(x, g=g, h=hvec[H])
  }
  
  df <- tibble(density = as.vector(ghpdf),
               x = rep(x, nh),
               h = factor(rep(hvec, each = nx)))
  
  out <- ggplot(df, aes(x = x, y = density, colour = h)) + theme_bw() +
    geom_line(linewidth = 1) +
    scale_colour_viridis_d(end = 0.9, option = "B") +
    scale_x_continuous(breaks = seq(-4, 5, 1)) +
    coord_cartesian(xlim = c(-4, 5)) +
    labs(x = "Observations",
         y = "Density") +
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = c(.1, .55),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(reverse=TRUE, override.aes = list(linewidth = 3)))
  
  out
}

