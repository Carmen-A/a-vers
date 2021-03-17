############################################################################################

                             # Anistro_Carmen_T3_P1.R #

############################################################################################

#Estructuras de selecci ??on
#Elabora un programa que, a partir de un vector de calificaciones de ex ??amenes parciales (3),
#y de un vector de calificaciones de tareas y ex ??amenes quincenales (10), determine si esa
#alumnx exentar ??ia o no el curso de bioinform ??atica con los criterios que establecimos 
#al inicio del curso que puedes consultar en el programa que en el classroom. Asume que en el 
#resto de los criterios tienen el 100%

#EP=50%
#EQ+T=30%
#PF=10
#Eart.=5%
#Part=5%

#exento > 70|7 = EP+EQ+T >|= 50|5

  EP<-c(8.5, 6.0, 4)
  EQT<-c(6, 9.75, 8, 4, 6.5, 7, 7, 8.5, 5, 9)
  
  x<-mean(EP)
  x
  y<-mean(EQT)
  y
  
  z<-50
  w<-30
  a<-100
  
  proporciónEP<-(x*z)/a
  proporción
  
  proporciónEQT<-(y*w)/a
  proporciónEQT
  
  t<-proporciónEP+proporciónEQT
  t
  
  if(t>5){print(paste("exento"))} else{print(paste("final"))}

  ############################################################################################
  
                               # Anistro_Carmen_T3_P2.R #
  
  ############################################################################################
  
  #Estructuras de selecci ??on
  #Elabora un programa que a partir de la fecha de nacimiento de cualquier persona, dando como 
  #datos de entrada d ??ia, mes y a ~no en formato num ??erico
  #calcule lo siguiente:
   #. Su edad en a ~nos
   #. Su signo zodiacal
   #. La estaci ??on del a ~no en el que naci ??o
  
  D<-22
  M<-10
  A<-1999
  
  #edad
  Aactual<- 2020
  E<-Aactual-A
  E

  #signo zodiacal
  #Acuario 20/1 - 18/2
  #Piscis 19/2 - 20/3
  #Aries 21/3 - 20/4
  #Tauro 21/4 - 20/5
  #Géminis 21/5 - 20/6
  #Cáncer 21/6 - 20/7
  #Leo 21/7 - 21/8
  #Virgo 22/8 - 22/9
  #Libra 22/9 - 22/10
  #Escorpio 23/10 - 22/11
  #Sagitario 23/11 - 20/12
  #Capricornio 21/12 - 19/1
  
  if((D>=20 && M==1)||(D<=18 && M==2)) 
  {print(paste("acuario"))} else if((D>=19 && M==2)||(D<=20 && M==3))
  {print(paste("piscis"))} else if((D>=21 && M==3)||(D<=20 && M==4)) 
  {print(paste("Aries"))} else if((D>=21 && M==4)||(D<=20 && M==5))
  {print(paste("tauro"))} else if((D>=21 && M==5)||(D<=20 && M==6)) 
  {print(paste("géminis"))} else if((D>=21 && M==6)||(D<=20 && M==7))
  {print(paste("cáncer"))} else if((D>=21 && M==7)||(D<=21 && M==8)) 
  {print(paste("leo"))} else if((D>=22 && M==8)||(D<=22 && M==9))
  {print(paste("virgo"))} else if((D>=22 && M==9)||(D<=22 && M==10)) 
  {print(paste("libra"))} else if((D>=23 && M==10)||(D<=22 && M==11))
  {print(paste("escorpio"))} else if((D>=23 && M==11)||(D<=20 && M==12)) 
  {print(paste("sagitario"))} else {print(paste("capricornio"))} 
   
  #estación del año
  #Primavera: 21/3 - 20/6
  #Verano: 21/6 - entre el 21/9
  #Otoño: 22/9 - 20/12.
  #Invierno: 21/12 - 20/3 
  D
  M
  if((D>=21 && M==3)||(D<=20 && M<=6)) {print(paste("primavera"))
    } else if((D>=21 && M==6)||(D<=21 && M<=9)) {print(paste("verano"))
    } else if((D>=22 && M==9)||(D<=20 && M<=12)) {print(paste("otoño"))
    } else {print(paste("invierno"))}

 ############################################################################################
  
                             # Anistro_Carmen_T3_P3.R #
  
  ############################################################################################
  
 # Alineamiento A partir de un archivo concatenado del genoma completo 10 virus elabora
  #un programa en R que realice los siguiente:
  
  library(Biostrings)
  cat<- readDNAStringSet("C:/Users/manis/Desktop/virus")
  
   # . Calcula el tama ~no de cada secuencia, es decir, el n ??umero de nucle ??otidos
  nchar(cat)
  
  # . Cambiar el nombre de las secuencias por uno m ??as corto, pro ejemplo, que incluya s ??olo
  #el nombre com ??un.
  names(cat)<- c("Abelson murine leukemia virus",  "Acanthocystis turfacea Chlorella",  
                 "Abaca bunchy top virus",  "Abalone shriveling syndrome-associated virus",
                 "Abutilon Brazil virus",  "Acanthamoeba polyphaga mimivirus", 
                 "Abutilon mosaic Bolivia virus",  "Abutilon mosaic Brazil", 
                 "Abalone herpesvirus", "Acartia tonsa copepod circovirus", "Abisko virus")
  head(cat)
  
  #. Generar un nuevo objeto para cada una de la siguientes operaciones
  #(a) El complemento
  catcome<-complement(cat)
  head(catcome)
  
  #(b) El reverso complemento
  CR<-reverseComplement(cat)
  CR
  reverse(catcome)
  
  #(c) El reverso
  R<-reverse(cat)
  R
  
  #(d) La secuencia traducida
  aa<-translate(cat)
  cat[[2]][288024] #hay un N en la posición 288024 de la seq 2
  
  IUPAC_CODE_MAP
  vmatchPattern("N", cat) #confirmadito
  chartr("N", "C", cat) # remplazar las N por C en el documento cat
  cat[[2]][288024] # checar :c whattt!

  aa<-translate(cat)# que todavía no
  
  aa<-translate(cat, no.init.codon =T, if.fuzzy.codon = "solve")
  aa
  
  #. Realiza un alineamiento por pares de las secuencias de amino ??acidos entre las dos
  #secuencias m ??as peque ~nas y las dos m ??as grandes.
  library(msa)
  length(aa[[1]])
  length(aa[[2]])
  length(aa[[3]])
  length(aa[[4]])
  length(aa[[5]])
  length(aa[[6]])
  length(aa[[7]])
  length(aa[[8]])
  length(aa[[9]])
  length(aa[[10]])
  length(aa[[11]])
  
  #secuencias más pequeñas:10 y 3
  #seucuencias más grandes: 6 y 2
  msa(aa, method = "Muscle")
  #. Encuentra los codones de inicio y de paro.
  
  
  #. Selecciona todas aquellas secuencias de nucle ??otidos que superen la media de sus lon-
   # gitudes y que mande a pantalla el nombre de c/u de ellas.
  
  
  
  
  