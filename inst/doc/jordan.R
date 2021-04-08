## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library("jordan")

## -----------------------------------------------------------------------------
x <- rrsm()  # "Random Real Symmetric Matrix"
y <- rrsm()  
z <- rrsm()  
x

## -----------------------------------------------------------------------------
x*100
x + y*3
x + 100

## -----------------------------------------------------------------------------
x*y

## -----------------------------------------------------------------------------
x*(y+z) - (x*y + x*z)

## -----------------------------------------------------------------------------
LHS <- x*(y*z)
RHS <- (x*y)*z
LHS-RHS

## -----------------------------------------------------------------------------
LHS <- (x*y)*(x*x)
RHS <- x*(y*(x*x))
LHS-RHS

## -----------------------------------------------------------------------------
M1 <- as.1matrix(x[1])
(M2 <- as.1matrix(x[2]))

## -----------------------------------------------------------------------------
(M1 %*% M2 + M2 %*% M1)/2 - as.1matrix(x[1]*x[2])

## -----------------------------------------------------------------------------
jj <- as.1matrix(x[1]*x[2])
jj-t(jj)

## -----------------------------------------------------------------------------
as.1matrix(rqhm(n=1,d=2))

## -----------------------------------------------------------------------------
x <- rqhm()
y <- rqhm()
(x*y)*(x*x) - x*(y*(x*x))

## -----------------------------------------------------------------------------
I <- rspin()
J <- rspin()
K <- rspin()
I
I*J - J*I   # commutative:
I*(J+K) - (I*J + I*K)  # distributive:
I*(J*K) - (I*J)*K  # not associative:
(I*J)*(I*I) - I*(J*(I*I))  # Obeys the Jordan identity

## -----------------------------------------------------------------------------
x <- ralbert()
y <- ralbert()
x
(x*y)*(x*x)-x*(y*(x*x)) # Jordan identity:

## ----label=define_U_and_diff--------------------------------------------------
U <- function(x){function(y){2*x*(x*y)-(x*x)*y}}
diff <- function(x,y,z){
     LHS <- U(x)(U(y)(U(x)(z)))
     RHS <- U(U(x)(y))(z)
     return(LHS-RHS)  # zero if Jacobson holds
}

## ----label=jacobsonverification,cache=TRUE------------------------------------
diff(ralbert(),ralbert(),ralbert())  # Albert algebra obeys Jacobson:
diff(rqhm(),rqhm(),rqhm()) # Quaternion Jordan algebra obeys Jacobson:
diff(rspin(),rspin(),rspin()) # spin factors obey Jacobson:

## ----label=defBH8G8-----------------------------------------------------------
B <- function(x,y,z){2*(x*(y*z) + (x*y)*z - (x*z)*y)}  # bracket function
H8 <- function(x,y,z){B(U(x)(U(y)(z)),z,x*y) - U(x)(U(y)(U(z)(x*y)))}
G8 <- function(x,y,z){H8(x,y,z)-H8(y,x,z)}

## ----label=verifyG8special,cache=TRUE-----------------------------------------
G8(rqhm(1),rqhm(1),rqhm(1))   # Quaternion Jordan algebra obeys G8:
G8(rspin(1),rspin(1),rspin(1)) # Spin factors obey G8:

## ----cache=TRUE---------------------------------------------------------------
G8(ralbert(1),ralbert(1),ralbert(1)) # Albert algebra does not obey G8:

## ----label=defineH9G9---------------------------------------------------------
L <- function(x){function(y){x*y}}
U <- function(x){function(y){2*x*(x*y)-(x*x)*y}}
U2 <- function(x,y){function(z){L(x)(L(y)(z)) + L(y)(L(x)(z)) - L(x*y)(z)}}
H9 <- function(x,y,z){2*U(x)(z)*U2(y,x)(U(z)(y*y)) - U(x)(U(z)(U2(x,y)(U(y)(z))))}
G9 <- function(x,y,z){H9(x,y,z)-H9(y,x,z)}

## ----verifyG9identity,cache=TRUE----------------------------------------------
G9(rqhm(1),rqhm(1),rqhm(1))  # Quaternion Jordan algebra obeys G9:

## ----albertH9G9,cache=TRUE----------------------------------------------------
G9(ralbert(1),ralbert(1),ralbert(1)) # Albert algebra does not obey G9:

