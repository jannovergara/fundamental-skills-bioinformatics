# title: "FSiB - Introduction to R"
# author: "Robert Lehmann" (modified "David Gomez-Cabrero")
# date: "2022"

### Week 1
### Data Types in R

## 0- Working in the right folder

getwd()
setwd()
 # you may want to replace this line 
 # with the outcome of the setwd() shown in your console

## 1- Assign scalar integer to variable
a <- 1
b <- 2

# Check class of variable a
class(a)

# perform operation on a and b
a*b
ab <- a*b

# Now to a character variable
c <- "Maria"
class(c)

# What happens when try to multiply an integer and a character?
a*c



## 2- Vectors

# Initialize vector directly.
v1 <- c(1,2,3,6,7,8)


# Initialize empty vector with type numeric
v.n <- vector('numeric', length = 3)
# Type can be any valid data type, numeric, logical, character, integer

# Can apply element-wise operations to vectors.
v1 + 1
v1 + 2
v1 * 2
v2 <- v1 * 2

# Obtain length of vector
length(v1)

# Obtain class of all elements in vector
class(v1)

# Element-wise operation of two vectors of equal length
v2 <- c(3:8)
v1 + v2

# Element-wise operation possible when length of one vector is a multiple of the other
v3 <- c(1,2,3)
v1 + v3 

# Error if not
A <- c(1,8,7,9)
B <- c(1,3,2)
A*B

# Vectors allow only elements of same data type
F <- c(1,"Gordon",3)

# Can add elements to existing vector
v2 <- c(v1, 15)
length(v1)
length(v2)



## 3- Lists

# Directly initialize list with varying data types
lst <- list(1, 2, 'three')

# Select a sublist of elements from list 
sb.lst <- lst[1]
class(sb.lst)

# Select a single element from list
elem <- lst[[1]]
class(elem)

# List elements can be named and accessed by name
nmd.lst <- list(first = 1, second = 2, third = 'three')
nmd.lst$second

nmd.lst <- list(x = 1, y = 2, z = 'three')
nmd.lst$y

# Get names of list elements
names(nmd.lst)


## 4- Matrices

# M1 <- matrix(fill, nrows, ncolums)
M1 <- matrix(0, 2, 3)

# Get dimensions of a matrix object
dim(M1)
nrow(M1)
ncol(M1)

# Inspect head of matrix
head(M1)
?head

M1
M1[1,1]

M1[1,1] <- 1
M1

# M2 and M3
M2 <- matrix(0,3,2)
M2[1,1] <- 3
M2[2,1] <- 2
M2[3,1] <- 6
M2[1,2] <- 1
M2[2,2] <- 3
M2[3,2] <- 2

M2 <- matrix( c(3,2,6,1,3,2) , 3 ,2)
M2

M3 <- matrix(c(2,2,3,6,1,1) ,2, 3) 
M3
M4 <- matrix(c(2,-1,3,-3,1,1) ,2, 3)

M2 + M3
M2 * M3
M2 + M4
M2 * M4
M4 + M3
M4 * M3

M4 * 3

##############
rownames(M4) <- c("patient 1","patient 2")
colnames(M4) <- c("He","Sugar","LDL")

M4[1,3]

M4[,3]
M4[,"LDL"]



## 5- Data Frame

v1 <- c(1,2,3)
v2 <- c(2,3,4)
v3 <- c(4,5,6)
v4 <- c("Gordon","Maria","Daniel")

df1 <- data.frame(col1 = v1, col2 = v2, col3 = v3, name = v4)

class(df1)
class(df1[,1])
class(df1[,4])

# Can access columns by name.
colnames(df1)

df1$name



## 6- Factor

# Initial data vector
data.vec <- c("small", "small", "medium", "large", "medium")
class(data.vec)

# Convert data to factor, all levels are considered equal, i.e. no order
data.factor <- factor(data.vec)
class(data.factor)
data.factor

# Can convert nominal factor to ordinal by introducing an order 
data.factor <- factor(data.vec, 
                      order = TRUE, 
                      levels =c("small", "medium", "large")) 
data.factor
?factor