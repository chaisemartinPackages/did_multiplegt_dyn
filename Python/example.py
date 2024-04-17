import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

import wooldridge

pandas2ri.activate()
wagepan = pandas2ri.py2ri(wooldridge.data('wagepan'))
#wagepan = pandas2ri.py2rpy(wooldridge.data('wagepan'))

DIDmultiplegtDYN = importr('DIDmultiplegtDYN')
did = DIDmultiplegtDYN.did_multiplegt_dyn(df = wagepan, outcome = 'lwage', group = 'nr', time = 'year', treatment = 'union', effects = 5.0, normalized = True, design = [1, "console"], controls = robjects.StrVector(["married", "hours"]))

# Display the results
print(did)

# Retrieve the ATE
print(did[1][3][0])