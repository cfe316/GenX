# Test Maintenance_Hours with maintenance based on core usage

This demonstrates the behavior of the fusion reactor maintenance constraint based on core usage.

The load is 
... 1
9   1
10  0
11  0
12  1
13  1
14  1
... 1

The reactor requries a 4-hour maintenance period every 20 hours of usage.
During the period of lower load, the core will go down for maintenance.

To try:
* Setting MAINT=0 will turn off the maintenance formulation. There should then be zero NSE.
* Turning off Integer_Commit for the core and changing Full_Power_Hours_Until_Maintenance => 86 should result in 0.25 cores down during the maintenance window. (This, rather than 92 hours, is the correct self-consistent solution.)

Misusage warnings
* Setting Full_Power_Hours_Until_Maintenance => 8 will result in two maintenance windows, but they will not be appropriately spaced. MAINT=2 is meant only for reactors that require maintenance less than once per simulation period.
