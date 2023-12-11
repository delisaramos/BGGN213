# class05
Delisa Ramos (PID:A69026881)

## Using GGPLOT :)

To use ggplot2, we first need to install it on our computers. To do this
we will is the function ‘install.packages()’. Before I use any package
functions I have to load them up with a ‘library()’ call, like so:

``` r
#install.packages("gifski")
#install.packages("gganimate")
#install.packages("gapminder")
#install.packages("patchwork")

library(gapminder)
library(ggplot2)
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ lubridate 1.9.3     ✔ tibble    3.2.1
    ✔ purrr     1.0.2     ✔ tidyr     1.3.0
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(gifski)
library(gganimate)
library(patchwork)
```

There is always the “base R” graphics system, i.e. ‘plot()’

``` r
plot(cars)
```

![](class05_files/figure-commonmark/cars%20in%20base%20R-1.png)

To use ggplot, I need to spell out at least 3 things: i. data (stuff I
want to plot) ii. aesthetics (aes) iii. geometries (geom\_)

``` r
ggplot(cars, aes(x=speed, y=dist)) + 
  geom_point() +
  geom_line()
```

![](class05_files/figure-commonmark/cars-1.png)

``` r
ggplot(cars, aes(x=speed, y=dist)) + 
  geom_point()
```

![](class05_files/figure-commonmark/cars-2.png)

``` r
#geom_smooth
ggplot(cars, aes(x=speed, y=dist)) + 
  geom_smooth(se=F, method = "lm") # lm = linear model
```

    `geom_smooth()` using formula = 'y ~ x'

![](class05_files/figure-commonmark/cars-3.png)

``` r
#add labs and color theme
#geom_smooth
ggplot(cars, aes(x=speed, y=dist)) + 
  geom_smooth() +
  ggtitle("Relationship between Speed and Distance") + 
  xlab("Speed") + ylab("Distance") +
  theme_bw()
```

    `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

![](class05_files/figure-commonmark/cars-4.png)

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

``` r
nrow(genes)
```

    [1] 5196

``` r
ncol(genes)
```

    [1] 4

``` r
colnames(genes)
```

    [1] "Gene"       "Condition1" "Condition2" "State"     

``` r
table(genes$State)
```


          down unchanging         up 
            72       4997        127 

``` r
round(127/5196, 2)
```

    [1] 0.02

``` r
g <- ggplot(genes, aes(x=Condition1, y= Condition2, col=State)) + 
  geom_point()
g + scale_color_manual(values=c("magenta", "navy", "purple")) + 
  ggtitle("Gene Expression Changes Upon Drug Treatment") +
  xlab("Control (no drugs)") + ylab("Drug Treatment")
```

![](class05_files/figure-commonmark/genes-1.png)

``` r
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"
gapminder <- read.delim(url)

head(gapminder)
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1952  28.801  8425333  779.4453
    2 Afghanistan      Asia 1957  30.332  9240934  820.8530
    3 Afghanistan      Asia 1962  31.997 10267083  853.1007
    4 Afghanistan      Asia 1967  34.020 11537966  836.1971
    5 Afghanistan      Asia 1972  36.088 13079460  739.9811
    6 Afghanistan      Asia 1977  38.438 14880372  786.1134

``` r
#filter to 2007
gapminder_07 <- gapminder %>% filter(year==2007)
colnames(gapminder_07)
```

    [1] "country"   "continent" "year"      "lifeExp"   "pop"       "gdpPercap"

``` r
ggplot(gapminder_07, aes(x=gdpPercap, y=lifeExp, col=continent, size=pop)) +
  geom_point(alpha=0.4)
```

![](class05_files/figure-commonmark/gapminder-1.png)

``` r
#alternatively
ggplot(gapminder_07, aes(x=gdpPercap, y=lifeExp, col=pop)) +
  geom_point(alpha=0.4)
```

![](class05_files/figure-commonmark/gapminder-2.png)

``` r
g07 <- ggplot(gapminder_07, aes(x=gdpPercap, y=lifeExp, col=continent, size=pop)) +
  geom_point(alpha=0.4) + 
  scale_size_area(max_size=10) + 
  ggtitle("Country's GDP percap vs. Life Expectancy", subtitle = "year: 2007")
g07
```

![](class05_files/figure-commonmark/gapminder-3.png)

``` r
# year 1957
gapminder_1957 <- gapminder %>% filter(year==1957)

g1957 <- 
ggplot(gapminder_1957, aes(x=gdpPercap, y=lifeExp, col=continent, size=pop)) +
  geom_point(alpha=0.7) + 
  scale_size_area(max_size=15) + 
  ggtitle("Country's GDP percap vs. Life Expectancy", subtitle = "year: 1957")
g1957
```

![](class05_files/figure-commonmark/gapminder-4.png)

``` r
both <- gapminder %>% filter(year==2007 | year==1957)
ggplot(both, aes(x=gdpPercap, y=lifeExp, col=continent, size=pop)) +
  geom_point(alpha=0.7) + 
  scale_size_area(max_size=10) + 
  ggtitle("Country's GDP percap vs. Life Expectancy", subtitle = "year: bith 1957 and 2007") + facet_wrap(~year)
```

![](class05_files/figure-commonmark/gapminder-5.png)

``` r
gapminder_top5 <- gapminder %>% 
  filter(year==2007) %>% 
  arrange(desc(pop)) %>% 
  top_n(5, pop)

gapminder_top5
```

            country continent year lifeExp        pop gdpPercap
    1         China      Asia 2007  72.961 1318683096  4959.115
    2         India      Asia 2007  64.698 1110396331  2452.210
    3 United States  Americas 2007  78.242  301139947 42951.653
    4     Indonesia      Asia 2007  70.650  223547000  3540.652
    5        Brazil  Americas 2007  72.390  190010647  9065.801

``` r
ggplot(gapminder_top5, aes(x=country, y=pop)) + 
  geom_col()
```

![](class05_files/figure-commonmark/bar%20charts-1.png)

``` r
colnames(gapminder_top5)
```

    [1] "country"   "continent" "year"      "lifeExp"   "pop"       "gdpPercap"

``` r
#life expectancy
ggplot(gapminder_top5, aes(x=country, y=lifeExp)) + 
  geom_col()
```

![](class05_files/figure-commonmark/bar%20charts-2.png)

``` r
#add color
ggplot(gapminder_top5, aes(x=country, y=lifeExp, fill=continent)) + 
  geom_col()
```

![](class05_files/figure-commonmark/bar%20charts-3.png)

``` r
ggplot(gapminder_top5, aes(x=country, y=lifeExp, fill=lifeExp)) + 
  geom_col()
```

![](class05_files/figure-commonmark/bar%20charts-4.png)

``` r
ggplot(gapminder_top5, aes(x=reorder(country, -pop), y=pop, fill=country)) + 
  geom_col()
```

![](class05_files/figure-commonmark/bar%20charts-5.png)

``` r
head(USArrests)
```

               Murder Assault UrbanPop Rape
    Alabama      13.2     236       58 21.2
    Alaska       10.0     263       48 44.5
    Arizona       8.1     294       80 31.0
    Arkansas      8.8     190       50 19.5
    California    9.0     276       91 40.6
    Colorado      7.9     204       78 38.7

``` r
#add new column... could also use mutate()
USArrests$State <- rownames(USArrests)

ggplot(USArrests) +
  aes(x=reorder(State,Murder), y=Murder) +
  coord_flip() + 
  geom_point() + #adds a point to end of each bar, helps with comparison
  geom_segment(aes(x=State, xend=State, 
                   y=0, yend= Murder), color="magenta")
```

![](class05_files/figure-commonmark/flipping%20bar%20charts-1.png)

``` r
# Setup nice regular ggplot of the gapminder data
ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  # Facet by continent
  facet_wrap(~continent) +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(year) +
  shadow_wake(wake_length = 0.1, alpha = FALSE)
```

``` r
#just messing around
(g07/g1957)
```

![](class05_files/figure-commonmark/facet%20wrap-1.png)
