## Info - 23 Apr 2024

This version of did_multiplegt_dyn is based off the 1.0.8 version available on CRAN. 
This version addresses rJava issues in some Linux and MacOS terminals by shutting down the the excel print features built in the **design** and **date_first_switch** options. 
These features were based on the 'xlsx' library, which depends on rJava. 
The output of these two options can be still saved in xlsx format using other libraries and retrieving the tables from the assigned object.