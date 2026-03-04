# Comtrade data

Eleven years of import and export data between 229 countries. The data
use the SITC Rev. 1 commodity classification, aggregated at the first
level (AG1).

## Format

A list consisting of a socioarray `Trade` and a vector `dollars2010` of
inflation rates. The socioarray gives yearly trade volume (exports and
imports) in dollars for 10 different commodity classes for eleven years
between 229 countries. This gives a five-way array. The first index is
the reporting country, so `Trade[i,j,t,k,1]` is what `i` reports for
exports to `j`, but in general this is not the same as
`Trade[j,i,t,k,2]`, what `j` reports as importing from `i`.

## Source

<https://comtrade.un.org/>, <https://www.measuringworth.com/>
