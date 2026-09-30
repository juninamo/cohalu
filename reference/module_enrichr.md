# Query enrichR for each co-localization module

\*\*Experimental.\*\* Convenience wrapper around \[enrichR::enrichr()\]
(needs the enrichR package and an internet connection): submits the gene
names of every module to the Enrichr web service and returns the
combined results. Only gene names are sent. For offline use, or to test
your own pathways or cell-type markers, use \[module_enrichment()\].

## Usage

``` r
module_enrichr(
  modules,
  databases = c("GO_Biological_Process_2023", "CellMarker_2024"),
  site = "Enrichr"
)
```

## Arguments

- modules:

  Output of \[colocalization_modules()\] or
  \[colocalization_submodules()\].

- databases:

  enrichR libraries, e.g. \`"GO_Biological_Process_2023"\`,
  \`"CellMarker_2024"\`, \`"MSigDB_Hallmark_2020"\`.

- site:

  Enrichr site passed to \[enrichR::setEnrichrSite()\].

## Value

A data.frame with a \`module\` and a \`database\` column added to the
enrichR output.
