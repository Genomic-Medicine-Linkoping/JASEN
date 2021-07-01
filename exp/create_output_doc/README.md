# Experimenting and testing output creation

- [Experimenting and testing output creation](#experimenting-and-testing-output-creation)
	- [Initialisation of `renv`](#initialisation-of-renv)
	- [Taking snapshot of the r environment](#taking-snapshot-of-the-r-environment)

This directory contains files related to experimenting and testing of final html report production. 

The testing environment was run inside conda environment defined by `env.yml` file and the r packages were installed inside a [`renv`](https://rstudio.github.io/renv/articles/renv.html)-environment. 

In order to get the environment running certain commands were necessary to have been run. They are described below.

## Initialisation of `renv`

```r
renv::init(bare = TRUE)
```

Due to an issue described [here](https://github.com/tidyverse/haven/issues/363#issuecomment-415331014), the following commands were necessary to have been run before installing `tidyverse` package:

```r
install.packages("withr")
withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("haven"), assignment = "+=")
withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("readxl"), assignment = "+=")
```
and installation of `tidyverse`:

```r
install.packages("tidyverse")
```

## Taking snapshot of the r environment

Finally, once everything was installed. A snapshot of the current state is in place.

By default, renv only writes out the packages you're currently using to the lockfile. If you want to just capture the state of the library without intersecting with the packages currently in use:

```r
renv::settings$snapshot.type("simple")
```

and the snapshot is taken with command:

```r
renv::snapshot()
```
