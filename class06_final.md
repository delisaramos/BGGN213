# class6
Delisa Ramos, pID A69026881

## Quarto

Quarto enables you to weave together content and executable code into a
finished document. To learn more about Quarto see <https://quarto.org>.

## ALl about functions in R

Every functinon in R has at least 3 things: - name - arguements (the
input(s) to your function) - the body

``` r
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ✔ purrr     1.0.2     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
student1 <- c(100, 100, 100, 100, 100, 100, 100, 90)
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
student3 <- c(90, NA, NA, NA, NA, NA, NA, NA)
```

Let’s start slow and find the average for student 1:

``` r
mean(student1)
```

    [1] 98.75

How can we drop the student’s lowest score?

``` r
min(student1)
```

    [1] 90

``` r
which.min(student1)
```

    [1] 8

``` r
student1 <- student1[-which.min(student1)]

#or

mean(student1[-which.min(student1)])
```

    [1] 100

will this work for student2?

``` r
mean(student2[-which.min(student2)], na.rm=T)
```

    [1] 92.83333

We can “mask” the NA or change them to be zero. THe rational here is if
you do not do a hw, then you get zero points.

``` r
student2
```

    [1] 100  NA  90  90  90  90  97  80

``` r
is.na(student2)
```

    [1] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE

``` r
student2[is.na(student2)] <- 0
```

``` r
grade <- function(student){
  #assigns any missing grades to 0
  student[is.na(student)] <- 0
  #takes the mean of the student's grades after removing the lowest score first
  mean(student[-which.min(student)])
}
```

We can now use this function called ‘grade’ to grade any student! yay

``` r
grade(student1)
```

    [1] 100

``` r
grade(student2)
```

    [1] 91

``` r
grade(student3)
```

    [1] 12.85714

``` r
student_hw <- read.csv("~/Desktop/student_homework.csv", row.names =1)

final_grades <- apply(student_hw, MARGIN=1, FUN=grade)

which.max(final_grades)
```

    student-18 
            18 

``` r
#18th student

#toughest hw

#assign all na assignments to 0
student_hw[is.na(student_hw)] <- 0
which.min(apply(student_hw, MARGIN = 2, FUN=mean))
```

    hw2 
      2 

``` r
which.min(apply(student_hw, MARGIN=2, FUN=sum))
```

    hw2 
      2 

``` r
#which assignment best correlates with overall final grade?
apply(student_hw, 2, cor, y=final_grades)
```

          hw1       hw2       hw3       hw4       hw5 
    0.4250204 0.1767780 0.3042561 0.3810884 0.6325982 

``` r
dim(mtcars)
```

    [1] 32 11

``` r
#install.packages("Bioconductor")
```
