# Experimenting and testing output creation


```r
renv::init(bare = TRUE)
```

By default, renv only writes out the packages you're currently using to the lockfile. If you want to just capture the state of the library without intersecting with the packages currently in use:

```r
renv::settings$snapshot.type("simple")
```

https://github.com/tidyverse/haven/issues/363#issuecomment-415331014

```r
install.packages("withr")
withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("haven"), assignment = "+=")
withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("readxl"), assignment = "+=")
```

```r
install.packages("tidyverse")
```
