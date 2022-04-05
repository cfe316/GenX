# Test Maintenance_Hours with Annual Cadence

This demonstrates the behavior of the fusion reactor maintenance constraint with an n-years cadence.
In the implementation, 'year' corresponds to the length of the time series, which here is a single day.

The load is 
... 4
9   4
10  3
11  3
12  4
13  4
14  4
... 4

The reactor requries maintenance one every 4 years, and the maintenance period is 4 hours.
There are 4 reactor cores, each of capacity 1 MW.
During the period of lower load, one of the four cores will go down for maintenance.

To try:
* Setting MAINT=0 will turn off the maintenance formulation. There should then be zero NSE.
* Setting Maintenance_Cadence_Years => 2 will result in two reactors needing to shut down at some point.
* Setting Maintenance_Cadence_Years => 1 will result in four reactors needing to shut down at some point.
