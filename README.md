# velocitree



### Fit models to data
```bash
# fit model to get estimated parameters and visualizations
velocitree fit \
	--tree test.nwk \
	--ridata crosses.csv \
	--spdata clades.csv \
	--models "all" \
	--functions "all" \
```


```
fitting 3 models (unpooled, pooled, partpooled)
fitting 5 functions (linear, logistic, quadratic, exponential, asymptotic)
```

### Perform model comparison
```bash

```


### Generate simulated data
```bash
# generate test datasets and visualizations
velocitree generate \
	--tree test.nwk \
	--model "all" \
	--function "all" \
```

```bash
# visualize 
```
