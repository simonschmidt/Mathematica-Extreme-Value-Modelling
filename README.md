Mathematica-Extreme-Value-Modelling
===================================
Under development

Example
-------

    Get["https://raw.github.com/simonschmidt/Mathematica-Extreme-Value-Modelling/master/ExtremeValues.m"]
    data = RandomVariate[MaxStableDistribution[30, 0.3, -0.3], 50];
    m = EVModelFit[data]
    m["Properties"]
    m["PlotDiagnostics"]

TODO 
----

* ::usage (and other documentation)
* Time-dependence
* Pareto (threshold)
* Friendly interface
* General touching up (Property renaming etc)

Done
----

* Figure out datastructure
* Probability/Quantile/ReturnLevel plots
* Confidence intervals based on normal approx
* Faster deviance statistic confidence interval

