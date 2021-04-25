sind<-function(d) sin(d*3.141593/180)
asind<-function(x) asin(x)/3.141593*180
cosd<-function(d) cos(d*3.141593/180)
acosd<-function(x) acos(x)*180/3.141593
tand<-function(d) tan(d*3.141593/180)
atand<-function(x) atan(x)/3.141593*180

rtd<-function(r) r/3.141593*180
dtr<-function(d) d/180*3.141593
#---------------------
Geom_wing<-function(x,...) UseMethod("Geom_wing")


Geom_AR<-function(x) UseMethod("Geom_AR")
Geom_AR.wing=function(x) x$ChordT - x$ChordR
Geom_Area = function(x) UseMethod("Geom_Area")
Geom_Area.wing=function(x){
  if(x$Type == 0|x$Type == 1) return(((x$ChordR+x$ChordT)*x$Span/2))
  #if(x$type == 2) #
}
Geom_lambda.wing = function(x) x$ChordT/x$ChordR
Geom_MAC.wing = function(x){
  l=Geom_lambda.wing(x)
  (1+l+l^2)/(1+L)
}
Geom_tanSwp.numeric<-function(x,n,m,tanSwpm,AR,l){
  tanSwpm-4/AR*(n-m)*(1-l)/(1+l)
}
Geom_tanSwp.wing<-function(x,n){
  Geom_tanSwp.numeric(x=x,n=n,m=0,tanSwpm = tand(x$Sweep),A = Geom_AR.wing(x),l=Geom_lambda.wing(x))
}
#-----------------------

#Stability and control simplified

SC_Vh=function(Sh,Lh,S,C) Sh*Lh/S/C

SC_B=function(lv,b,Gamma,Cl) lv/b*Gamma/Cl

SC_simp<-function(x,...) UseMethod("SC_simp")
SC_simp.conventionalConcept=function(x,Xcg,Xnp = NULL,gamma.include =F,Cl = NA,Gamma=NA){
  MACM = MAC(x$WM);MACH= MAC(x$WH);MACV=MAC(x$WV)
  #absolute position of m,h,v
  Pos = list(xM = MACM$AF+x$WM$x,
             xH = MACH$AF+x$WH$x,
             xV = MACV$AF+x$WV$x)
  Lever = list(L = Xcg - Pos$xM,
               LH= Xcg - Pos$xH,
               LV= Xcg - Pos$xV)
  Vh = Geom_Area(x$WH)*Lever$LH / Geom_Area(x$WM) / MACM$ChordAvg#Sh Lh/S/c
  Vv = Geom_Area(x$WV)*Lever$LV / Geom_Area(x$WV) / MACM$ChordAvg
  if (is.null(Xnp)){
    Xnp_OCW_est<-MACM$ChordAvg*(.25+(1+2/Geom_AR(x$WM))/(1+2/Geom_AR(x$WH)))*(1-4/(2+Geom_AR(x$WM))*Vh)+x$WM$x
    #message(paste0('No Xnp Detected, using estimation from OCW Lab8 Notes P4',Xnp_OCW_est))
    Xnp = Xnp_OCW_est
  }
  SM =(Xcg - Xnp)/ MACM$ChordAvg#;message('OK')#Static Margin
  #Nondimensionization (torque d(T)/d(alpha) / (.5*Cl*rho*v^2) )
  M_de<-Lever$L  / MACM$ChordAvg + Vh#l/b - St*L1/(Sw*b) P67 fuyangwendingxishu
  Coef = cbind(Vh,Vv,SM,M_de)
  
  if(gamma.include==T){
    B = Lever$LV / x$WM$Span * Gamma / Cl#Gamma in Deg Blaine Rowdon Factor of spiral stability
    VvB = Vv*B
    Coef = cbind(Coef,B)
    Coef = cbind(Coef,VvB)
  }
  
  Raw = list(Xcg=Xcg,Xnp=Xnp)
  re = list(From =deparse(substitute(x)),Raw = Raw,Pos = Pos,Lever = Lever,Coef=Coef)
  class(re) = 'SCSimpOut'
  re
}

print.SCSimpOut<-function(x){
  print(paste('Simplified Stability and Control Analysis for',x$From))
  print(paste('Center of Gravity (Xcg):',x$Raw$Xcg,'|Neutral Point(Xnp):',x$Raw$Xnp))
  print('...Coefficients:')
  print(x$Coef)
}

#---------------------------------------
#Naming refers to\

#  Root Chord (A):
#MAC Graphic
#  Tip Chord (B):
#  Sweep Distance (S):
#  Half Span (Y):
#  C  of Gravity (CG):
#Sweep Distance @ MAC (C) =
#Mean Aerodynamic Chord (MAC) =
#MAC Distance from Root (d) =
#http://www.nasascale.org/p2/wp-content/uploads/mac-calculator.htm
#
#MAC
MAC_S<-function(span,sweep) span * tand(sweep)/2
MAC_C<-function(S,A,B) (S*(A+2*B)) / (3*(A+B))#
MAC_MAC<-function(A,B) A-(2*(A-B)*(0.5*A+B) / (3*(A+B)))
MAC_d<-function(A,B,Y)  (2*Y*(0.5*A+B)) / (3*(A+B))
#MAC_B.P.<-function(CG,C,MAC) ((CG-C) / MAC)
MAC_AF<-function(B.P.,MAC,C) C + B.P.* MAC

MAC<-function(x,...) UseMethod('MAC')


MAC.wing<-function(wing,B.P=NA){
  if(is.na(B.P)){
    try(B.P<-wing$Xf_C)
    message('B.P is missing from call, Extracting from wing')
  }
  if(wing$Type==0){
    re<-list(
      xSweep = MAC_S(wing$Span,wing$Sweep),#Sweep Distance @ tip
      xMAC  = MAC_C(MAC_S(wing$Span,wing$Sweep),wing$ChordR,wing$ChordT),#Sweep Distance @ MAC
      ChordAvg = MAC_MAC(wing$ChordR,wing$ChordT),
      yMAC = MAC_d(wing$ChordR,wing$ChordT,wing$Span/2),
      AF=MAC_AF(B.P,MAC_MAC(wing$ChordR,wing$ChordT), MAC_C(MAC_S(wing$Span,wing$Sweep),wing$ChordR,wing$ChordT)))
    re$Type=0
    class(re)<-c('MAC','list')
    re
  }else if(wing$Type==1){
    wing$Span<-wing$Span*2
    wing$Type<-0
    #wing$
    rr=MAC(wing)
    rr$Type=1
    return(rr)
  }else if(wing$Type==2){
    
  }
}


#MAC_CG<-function(B.P.,MAC,C) C + B.P.* MAC
#afprop<-MAC(wing.default,B.P = .25)


pointsxy.MAC<-function(afprop,xf=0){
  lines(c(xf+afprop$xMAC,xf+afprop$xMAC+afprop$ChordAvg),c(afprop$yMAC,afprop$yMAC),type='l',col='red')
  points(xf+afprop$AF,afprop$yMAC,col='red')
}

#w = list(x=0,ChordR = .5,ChordT = .3,Sweep = 12,Span = 1.5)
#class(w) = 'wing'

#MAC.wing(w,Xcg = .21)
#
plotxy<-function(x,...) UseMethod('plotxy')
pointsxy<-function(x,...) UseMethod('pointsxy')


plotxy.concept<-function(concept,asp = 1){
  plot(c(0,FU$Length),c(0,FU$Length),type = 'n',asp = asp)
  plotxy.canard<-function(concept){
    pointsxy(WM)
    pointsxy(WC)
  }
  plotxy.canard2v<-function(concept){
    pointsxy(WM)
    pointsxy(WC)
  }
  plotxy.conventional<-function(concept){
    pointsxy(WM)
    pointsxy(WH)
  }
  plotxy.conventional2v<-function(concept){
    pointsxy(WM)
    pointsxy(WH)
  }
  plotxy.aerobatics<-function(concept){
    pointsxy(WM)
    pointsxy(WH)
  }
  plotxy(concept)
}

pointsxy.concept<-function(concept){
  plotxy.canard<-function(concept){
    pointsxy(WM)
    pointsxy(WC)
  }
  plotxy.canard2v<-function(concept){
    pointsxy(WM)
    pointsxy(WC)
  }
  plotxy.conventional<-function(concept){
    pointsxy(WM)
    pointsxy(WH)
  }
  plotxy.conventional2v<-function(concept){
    pointsxy(WM)
    pointsxy(WH)
  }
  plotxy.aerobatics<-function(concept){
    pointsxy(WM)
    pointsxy(WH)
  }
  pointsxy(concept)
}

#pointsxy.MAC<-function(MAC,x=0){
#   points(c(x + MAC$D_Sweep_atMAC,x+ MAC$D$Sweep_atMAC +MAC$ChordAvg),c(MAC$yMAC,MAC$yMAC),type = 'l',col = 'red')
#}

extend<-function(x,...) UseMethod('extend')

#------------------
plotxy<-function(x,...) UseMethod('plotxy')

plotxy.wing<-function(wing,add=F,col='black',lty=1){
  if(add==T){
    pointsxy.wing(wing,col=col,lty=lty)
  }else{
    #if(missing(wing$Dihedral)) wing$Dihedral=0
    span<-wing$Span/2
    x=wing$x;y=wing$y;z=wing$z
    front<-c(x,y,z)
    #dim(front) = c(1,3)
    rear<-c(x+wing$ChordR,y,z*wing$Dihedral)
    #dim(rear) = c(1,3)
    frontT<-c(x+span*sind(wing$Sweep),y+span*cosd(wing$Sweep),z+span*cosd(wing$Dihedral))
    #dim(frontT) = c(1,3)
    rearT<-c(x+span*sind(wing$Sweep)+wing$ChordT,y+span*cosd(wing$Sweep),z+span*cosd(wing$Dihedral))
    #dim(rearT) = c(1,3)
    pointCloud<-rbind(as.numeric(front)[1:2],as.numeric(frontT)[1:2],as.numeric(rearT)[1:2],as.numeric(rear)[1:2])
    plot(pointCloud[,1],pointCloud[,2],xlab='x',ylab='y',type='l',asp=1,col=col,lty=lty)
  }
}

pointsxy.wing<-function(wing,col='black',lty=1){
  #if(missing(wing$Dihedral)) wing$Dihedral=0
  span<-wing$Span/2
  x=wing$x;y=wing$y;z=wing$z
  front<-c(x,y,z)
  rear<-c(x+wing$ChordR,y,z*wing$Dihedral)
  frontT<-c(x+span*sind(wing$Sweep),y+span*cosd(wing$Sweep),z+span*cosd(wing$Dihedral))
  rearT<-c(x+span*sind(wing$Sweep)+wing$ChordT,y+span*cosd(wing$Sweep),z+span*cosd(wing$Dihedral))
  pointCloud<-rbind(as.numeric(front)[1:2],as.numeric(frontT)[1:2],as.numeric(rearT)[1:2],as.numeric(rear)[1:2])
  points(pointCloud[,1],pointCloud[,2],xlab='x',ylab='y',type='l',asp=1,col=col,lty=lty)
}

plotxy.fuse<-function(fuse,add=F,col='black',lty=1){
  
  if(add==T){
    points(fuse$x,fuse$r,col=col,type='l',lty=lty)
  }else{
    plot(fuse$x,fuse$r,xlab='x',ylab='y',type='l',asp=1,col=col,lty=lty)
  }
}
pointsxy.fuse<-function(fuse,col='black',lty=1){
  points(fuse$x,fuse$r,col=col,type='l',lty=lty)
}

#plotxy(fuselage.default,add=T)
plotxy.conventionalConcept = function(concept){
  #range of plotting
  plot(c(0,max(concept$fuselage$x)),c(0,max(concept$WM$Span/2)),type='n',asp=1,
       xlab = 'x',ylab = 'y')
  pointsxy(concept$fuselage)
  pointsxy(concept$WM)
  pointsxy(concept$WH)
  pointsxy(concept$WV,col = 'blue',lty=2)
  try(pointsxy(MAC(concept$WM),xf=concept$WM$x))
  try(pointsxy(MAC(concept$WH),xf=concept$WH$x))
  try(pointsxy(MAC(concept$WV),xf=concept$WV$x))
}






#-----------------------------------


library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Conventional Concept simplified"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      #WM
      sliderInput("Xcg",
                  "Xcg[cm]",
                  min=1,
                  max=120,
                  value = 25),
      sliderInput("WM.Span",
                  "WM.Span[cm]",
                  min = 1,
                  max = 200,
                  value =120),
      sliderInput("WM.Sweep",
                  "WM.Sweep[deg]",
                  min = 1,
                  max = 40,
                  value = 10),
      sliderInput("WM.ChordR",
                  "WM.ChordR[cm]",
                  min = 18,
                  max = 60,
                  value = 20),
      sliderInput("WM.ChordT",
                  "WM.ChordT[cm]",
                  min = 14,
                  max = 60,
                  value = 10),
      sliderInput("WM.Xf_C",
                  "WM.Xf_C[%]",
                  min = 1,
                  max = 100,
                  value = 25),
      sliderInput("WM.x",
                  "WM.x[cm]",
                  min = 1,
                  max = 200,
                  value = 20),
      sliderInput("WM.y",
                  "WM.y[cm]",
                  min = -100,
                  max = 100,
                  value = 0),
      sliderInput("WM.z",
                  "WM.z[cm]",
                  min = -100,
                  max = -100,
                  value = 0),#WH
      sliderInput("WH.Span",
                  "WH.Span[cm]",
                  min = 1,
                  max = 200,
                  value = 15),
      sliderInput("WH.Sweep",
                  "WH.Sweep[deg]",
                  min = 1,
                  max = 40,
                  value = 10),
      sliderInput("WH.ChordR",
                  "WH.ChordR[cm]",
                  min = 1,
                  max = 60,
                  value = 10),
      sliderInput("WH.ChordT",
                  "WH.ChordT[cm]",
                  min = 1,
                  max = 60,
                  value = 5),
      sliderInput("WH.Xf_C",
                  "WH.Xf_C[%]",
                  min = 1,
                  max = 100,
                  value = 25),
      sliderInput("WH.x",
                  "WH.x[cm]",
                  min = 1,
                  max = 200,
                  value = 80),
      sliderInput("WH.y",
                  "WH.y[cm]",
                  min = -100,
                  max = 100,
                  value = 0),
      sliderInput("WH.z",
                  "WH.z[cm]",
                  min = -100,
                  max = 100,
                  value = 0),#WV
      sliderInput("WV.Span",
                  "WV.Span[cm]",
                  min = 1,
                  max = 200,
                  value = 15),
      sliderInput("WV.Sweep",
                  "WV.Sweep[deg]",
                  min = 1,
                  max = 40,
                  value = 10),
      sliderInput("WV.ChordR",
                  "WV.ChordR[cm]",
                  min = 1,
                  max = 60,
                  value = 10),
      sliderInput("WV.ChordT",
                  "WV.ChordT[cm]",
                  min = 1,
                  max = 60,
                  value = 5),
      sliderInput("WV.Xf_C",
                  "WV.Xf_C[%]",
                  min = 1,
                  max = 100,
                  value = 25),
      sliderInput("WV.x",
                  "WV.x[cm]",
                  min = 1,
                  max = 200,
                  value = 80),
      sliderInput("WV.y",
                  "WV.y[cm]",
                  min = -100,
                  max = 100,
                  value = 0),
      sliderInput("WV.z",
                  "WV.z[cm]",
                  min = -100,
                  max = -100,
                  value = 0),
      textInput("Length","Length",value = '1.05'),
      textInput("r","r","c( 0.050,0.075,0.075,0.020,0.010,0.000)"),
      textInput("x","x","c(0.00,0.20,0.40,0.90,1.00,1.05)"),
      width =3
    ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot"),
      tabsetPanel(type = "tabs",
                  tabPanel("Main Wing MAC", tableOutput("wm")),
                  tabPanel("Horz tail MAC", tableOutput('wh')),
                  tabPanel("Vert tail MAC", tableOutput('wv'))),
      tableOutput("SC")
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$wm = renderTable({
    WM = list(Span = input$WM.Span/100,
              Sweep = input$WM.Sweep,
              ChordR = input$WM.ChordR/100,
              ChordT = input$WM.ChordT/100,
              Type = 0,
              Xf_C = input$WM.Xf_C/100,
              x = input$WM.x/100,
              y = input$WM.y/100,
              z = input$WM.z/100)
    class(WM) = 'wing'
    MAC(WM)
  })
  output$wh = renderTable({
    WH = list(Span = input$WH.Span/100,
              Sweep = input$WH.Sweep,
              ChordR = input$WH.ChordR/100,
              ChordT = input$WH.ChordT/100,
              Type = 0,
              Xf_C = input$WH.Xf_C/100,
              x = input$WH.x/100,
              y = input$WH.y/100,
              z = input$WH.z/100)
    class(WH) = 'wing'
    MAC(WH)
  })
  output$wv = renderTable({
    WV = list(Span = input$WV.Span/100,
              Sweep = input$WV.Sweep,
              ChordR = input$WV.ChordR/100,
              ChordT = input$WV.ChordT/100,
              Type = 0,
              Xf_C = input$WV.Xf_C/100,
              x = input$WV.x/100,
              y = input$WV.y/100,
              z = input$WV.z/100)
    class(WV) = 'wing'
    MAC(WV)
  })
  output$SC = renderTable({
    
    WM = list(Span = input$WM.Span/100,
              Sweep = input$WM.Sweep,
              ChordR = input$WM.ChordR/100,
              ChordT = input$WM.ChordT/100,
              Type = 0,
              Xf_C = input$WM.Xf_C/100,
              x = input$WM.x/100,
              y = input$WM.y/100,
              z = input$WM.z/100)
    class(WM) = 'wing'
    WH = list(Span = input$WH.Span/100,
              Sweep = input$WH.Sweep,
              ChordR = input$WH.ChordR/100,
              ChordT = input$WH.ChordT/100,
              Type = 0,
              Xf_C = input$WH.Xf_C/100,
              x = input$WH.x/100,
              y = input$WH.y/100,
              z = input$WH.z/100)
    class(WH) = 'wing'
    WV = list(Span = input$WV.Span/100,
              Sweep = input$WV.Sweep,
              ChordR = input$WV.ChordR/100,
              ChordT = input$WV.ChordT/100,
              Type = 1,
              Xf_C = input$WV.Xf_C/100,
              x = input$WV.x/100,
              y = input$WV.y/100,
              z = input$WV.z/100)
    class(WV) = 'wing'
    
    fuse = list(Length = as.numeric(input$Length),
                x = eval(parse(text = (input$x))),
                r = eval(parse(text = (input$r))))
    class(fuse) = "fuse"
    x = list(WM=WM,WH=WH,WV=WV,fuselage=fuse)
    class(x$WM) = 'wing'
    class(x$WH) = 'wing'
    class(x$WV) = 'wing'
    class(x) = "conventionalConcept"
    #plotxy(plane_init)
    s = SC_simp(x,input$Xcg/100)
    re = as.vector(c(unlist(s$From),unlist(s$Raw),unlist(s$Pos), unlist(s$Lever),unlist(s$Coef)))
    
    dim(re) = c(1,13)
    colnames(re) = c("From","Xcg","Xnp","Pos.xM","Pos.xH","Pos,xV",
                     "Lever","Lever.H","Lever.V",
                     "Coef.Vh","Coef.Vv","SM","M_de")
    re
  })
  output$distPlot <- renderPlot({
    WM = list(Span = input$WM.Span/100,
              Sweep = input$WM.Sweep,
              ChordR = input$WM.ChordR/100,
              ChordT = input$WM.ChordT/100,
              Type = 0,
              Xf_C = input$WM.Xf_C/100,
              x = input$WM.x/100,
              y = input$WM.y/100,
              z = input$WM.z/100)
    class(WM) = 'wing'
    WH = list(Span = input$WH.Span/100,
              Sweep = input$WH.Sweep,
              ChordR = input$WH.ChordR/100,
              ChordT = input$WH.ChordT/100,
              Type = 0,
              Xf_C = input$WH.Xf_C/100,
              x = input$WH.x/100,
              y = input$WH.y/100,
              z = input$WH.z/100)
    class(WH) = 'wing'
    WV = list(Span = input$WV.Span/100,
              Sweep = input$WV.Sweep,
              ChordR = input$WV.ChordR/100,
              ChordT = input$WV.ChordT/100,
              Type = 1,
              Xf_C = input$WV.Xf_C/100,
              x = input$WV.x/100,
              y = input$WV.y/100,
              z = input$WV.z/100)
    class(WV) = 'wing'
    
    fuse = list(Length = as.numeric(input$Length),
                x = eval(parse(text = (input$x))),
                r = eval(parse(text = (input$r))))
    class(fuse) = "fuse"
    x = list(WM=WM,WH=WH,WV=WV,fuselage=fuse)
    class(x$WM) = 'wing'
    class(x$WH) = 'wing'
    class(x$WV) = 'wing'
    class(x) = "conventionalConcept"
    #plotxy(plane_init)
    plotxy.conventionalConcept(x)#
    abline(v = input$Xcg/100,col = 'red')
    
    s = SC_simp(x,input$Xcg/100)
    re = as.vector(c(unlist(s$From),unlist(s$Raw),unlist(s$Pos), unlist(s$Lever),unlist(s$Coef)))
    
    dim(re) = c(1,13)
    colnames(re) = c("From","Xcg","Xnp","Pos.xM","Pos.xH","Pos,xV",
                     "Lever","Lever.H","Lever.V",
                     "Coef.Vh","Coef.Vv","SM","M_de")
    abline(v = re[1,3],col = "blue")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

