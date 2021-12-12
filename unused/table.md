
|true | est|   n|      freq|
|:----|---:|---:|---------:|
|BRCA |   1| 148| 0.1847690|
|BRCA |   2| 152| 0.1897628|
|COAD |   4|  78| 0.0973783|
|KIRC |   3|   2| 0.0024969|
|KIRC |   4|   3| 0.0037453|
|KIRC |   5| 141| 0.1760300|
|LUAD |   1|   2| 0.0024969|
|LUAD |   3| 107| 0.1335830|
|LUAD |   4|  32| 0.0399501|
|PRAD |   1|   1| 0.0012484|
|PRAD |   3| 135| 0.1685393|

### 1. spcaRcpp (alpha = beta = 1e-4) = 0.897628
### THIS IS IT
|true | est|   n|      freq|
|:----|---:|---:|---------:|
|BRCA |   2|   1| 0.0012484|
|BRCA |   3| 299| 0.3732834|
|COAD |   2|   3| 0.0037453|
|COAD |   4|  75| 0.0936330|
|KIRC |   2|   1| 0.0012484|
|KIRC |   5| 145| 0.1810237|
|LUAD |   2| 139| 0.1735331|
|LUAD |   3|   2| 0.0024969|
|PRAD |   1| 134| 0.1672909|
|PRAD |   3|   2| 0.0024969|

### 2. spcaRcpp (alpha = beta = 0) = 0.8938827
|true | est|   n|      freq|
|:----|---:|---:|---------:|
|BRCA |   1| 207| 0.2584270|
|BRCA |   2|  93| 0.1161049|
|COAD |   4|  78| 0.0973783|
|KIRC |   4|   1| 0.0012484|
|KIRC |   5| 145| 0.1810237|
|LUAD |   2|   2| 0.0024969|
|LUAD |   4| 139| 0.1735331|
|PRAD |   3| 132| 0.1647940|
|PRAD |   4|   4| 0.0049938|


### spca + kmeans_cluster
|true | est|   n|      freq|
|:----|---:|---:|---------:|
|BRCA |   1| 250| 0.3121099|
|BRCA |   2|  50| 0.0624220|
|COAD |   2|  78| 0.0973783|
|KIRC |   2|  90| 0.1123596|
|KIRC |   5|  56| 0.0699126|
|LUAD |   2|   2| 0.0024969|
|LUAD |   4| 139| 0.1735331|
|PRAD |   3| 136| 0.1697878|


### spca + kmeans_cluster (method = "gkmeans++")
|true | est|   n|      freq|
|:----|---:|---:|---------:|
|BRCA |   2| 300| 0.3745318|
|COAD |   2|   2| 0.0024969|
|COAD |   4|  76| 0.0948814|
|KIRC |   2|   1| 0.0012484|
|KIRC |   5| 145| 0.1810237|
|LUAD |   1| 136| 0.1697878|
|LUAD |   2|   5| 0.0062422|
|PRAD |   3| 136| 0.1697878|