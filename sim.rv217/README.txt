
Thu Sep 17 17:02:09 PDT 2015
Chris Warth

Simulated HIV samples mimicking the RV217 samples distributed as part
of the play data for the Gates HIV-1 founder inference bakeoff.

The idea here was to simulate HIV-1 samples that have only about
10-12 samples per timepoint and see how various methods perform when
inferring founder sequences.

The SConstruct file is the heart of the process for generating the
simulted sequences.  Scripts in the analysis/ subdirectory are used to
measure the distance from inferred to actual founder and plot the
results.
