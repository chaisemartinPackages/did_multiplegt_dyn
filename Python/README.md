# Running did_multiplegt_dyn in Python

The **DIDmultiplegtDYN** library has been developed for the R language. As a result, one could use the **rpy2** library to load the package in a Python environment. All the most updated options of the did_multiplegt_dyn command can be accessed, even though some tweaks are needed to adjust Python classes such that input objects can be read correcly. Lastly, the output can be assigned to a Python object and parsed with nested list syntax.

## The rpy2 library
The **rpy2** enables R commands and packages in Python scripts. It can be installed by running
```r
pip install rpy2
```
For Windows users, it is best to install rpy2 using **wheel**. First, get wheel:
```r
pip install wheel
```
Then, get the latest rpy2 .whl extension [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/#rpy2). Download the .whl file and open **cmd** in the directory where the .whl file is located. Run:
```r
pip install `whl_file'.whl
```
replacing `whl_file' with the name of the file. Now you can check that rpy2 is installed by loading python and checking whether **import rpy2** returns errors.

Before using rpy2 in Windows, you should also update your R PATH information. Browse *edit environmental variables* in the search menu and do the following:
1. Add R to PATH (with R 4.3.3, you should add a bin path similar to "C:\Program Files\R\R-4.3.3\bin\x64")
2. Create a new variable called R_HOME and assign it to the R version directory (e.g., "C:\Program Files\R\R-4.3.3")
3. Create a new variable called R_USER and assign it to the entry in the *user* field from running Sys.info() in R.

## R-Python integration
Start by loading the required libraries:
```r
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
```
We will need *robjects* to load string vectors, *importr* to load DIDmultiplegtDYN and *pandas2ri* to convert a pandas dataframe into an R dataframe. The last one should be explicitely activated in the script. To showcase the routine, we load the dataset from Vella and Verkeek ([1998](https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1099-1255(199803/04)13:2%3C163::AID-JAE460%3E3.0.CO;2-Y)):
```r
import wooldridge

pandas2ri.activate()
wagepan = pandas2ri.py2ri(wooldridge.data('wagepan'))
```

If you installed rpy2 using just pip install rpy2, then replace the last line with
```r
wagepan = pandas2ri.py2rpy(wooldridge.data('wagepan'))
```
Once the dataframe is loaded, we can set up the estimation options. The arguments of the options have different classes (numeric, list, character, ...). These R classes are normally different than Python ones ("in R, *almost* everything is a vector"), but a bit of tweaking is enough to ensure that all the options of **did_multiplegt_dyn** can be correctly assigned. Consider the following equivalences:
+ string arguments ("outcome", "group", "time", "treatment", "by", "cluster", "weight", "switchers", "save_results") can be specified as Python strings (e.g. outcome = 'lwage');
+ integer arguments ("effects","placebo", "ci_level", "continuous", "bootstrap", "by_path") can be specified as floats with zero decimal part (e.g. effects = 3.0) - This is due to the fact that the R integer class is quite restrictive, so the syntax check looks for numeric objects such that the remainder from 1 is 0;
+ list arguments ("predict_het", "design", "date_first_switch") can be specified as Python lists (e.g. design = [1, "console"]);
+ string vector arguments (e.g., arguments that can be strings or vectors of strings, that is, "controls", "trends_nonparams") can be specified as a single string in the single variable case or using *StrVector* in case of multiple variables (e.g. controls = robjects.StrVector(["married", "black"]));
+ logical arguments ("normalized", "normalized_weights", "effects_equal", "trends_lin", "same_switchers", "same_switchers_pl", "graph_off", "save_sample", "less_conservative_se", "dont_drop_larger_lower", "drop_if_d_miss_before_first_switch") can be specified using True/False.

The results can be assigned to an object that has it own print method. The base syntax for this object is the same as R lists, which implies that the object can be accessed through subsetting.

## Example
```r
DIDmultiplegtDYN = importr('DIDmultiplegtDYN')
did = DIDmultiplegtDYN.did_multiplegt_dyn(df = wagepan, outcome = 'lwage', group = 'nr', time = 'year', treatment = 'union', effects = 5.0, normalized = True, design = [1, "console"], controls = robjects.StrVector(["married", "hours"]))

# Display the results
print(did)

# Retrieve the ATE
print(did[1][3][0])
```

