[![Build Status](https://travis-ci.org/trosendal/paratb.svg?branch=master)](https://travis-ci.org/trosendal/paratb)
# paratb

This is an R package of an infectious disease model of
paratuberculosis in cattle.

## Installation

Use the `remotes` package to install the package and its dependency from
github like this:

```{r}
library("remotes")
install_github("trosendal/paratb")
```

## Usage

The usage of the model is described in a vignette that accompanies the
package. After installing the package, you should be able to view the
vignette web page like this:

```{r}
library(paratb)
vignette("model", package = "paratb")
```
