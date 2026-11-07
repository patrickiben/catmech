# Read a strict YAML assumption specification

Read a strict YAML assumption specification

## Usage

``` r
read_catmech_spec(path)
```

## Arguments

- path:

  Path to a YAML file. YAML expression evaluation is disabled.

## Value

A validated catmech_spec object. Legacy prototype YAML is deliberately
rejected; see inst/specs/lab.yaml for the supported schema.

## Examples

``` r
spec_file <- system.file("specs", "lab.yaml", package = "catmech")
if (file.exists(spec_file)) {
  sp <- read_catmech_spec(spec_file)
  print(sp)
}
#> $version
#> [1] "1"
#> 
#> $value
#> [1] "lab_value"
#> 
#> $category
#> [1] "lab_category"
#> 
#> $levels
#> [1] "Low"    "Medium" "High"  
#> 
#> $thresholds
#> [1] -Inf  5.5  7.0  Inf
#> 
#> $flag
#> [1] "outlier_flag"
#> 
#> $flag_yes
#> [1] "Yes"
#> 
#> $flag_no
#> [1] "No"
#> 
#> $corruption_probability
#> [1] 1
#> 
#> $conflict_policy
#> [1] "reconstruct"
#> 
#> $min_donors
#> [1] 1
#> 
#> $n_donors
#> [1] 5
#> 
#> $fallback
#> [1] "error"
#> 
#> attr(,"class")
#> [1] "catmech_spec"
```
