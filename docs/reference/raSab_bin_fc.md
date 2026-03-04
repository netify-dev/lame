# Simulate a and Sab from full conditional distributions under bin likelihood

Simulate a and Sab from full conditional distributions under bin
likelihood

## Usage

``` r
raSab_bin_fc(Z, Y, a, b, Sab, Sab0=NULL, eta0=NULL, SS = round(sqrt(nrow(Z))))
```

## Arguments

- Z:

  a square matrix, the current value of Z

- Y:

  square binary relational matrix

- a:

  current value of row effects

- b:

  current value of column effects

- Sab:

  current value of Cov(a,b)

- Sab0:

  prior (inverse) scale matrix for the prior distribution

- eta0:

  prior degrees of freedom for the prior distribution

- SS:

  number of iterations

## Value

- Z:

  new value of Z

- Sab:

  new value of Sab

- a:

  new value of a

## Author

Peter Hoff
