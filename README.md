# Monography
Scrips that i have used in my monography


==================================================================
CMS OPEN DATA 2011
==================================================================

JpsiAnalyzerOpen2011.cc
It is a file that work on CMSSW. This file get the muons cinematc variables with a trigger and put in a tree.

JpsiAnalyzerOpen2011_2.cc
It is the almost the same file above, but this one get the muons cinematc variables whithou trigger. They have diffent name because the GitHub does not keep file with same name.

JPSIAnaOpenData2011_Model_Trigger.py
This file work with "JpsiAnalyzerOpen2011.cc". This script in python was made to get acess the online root files in CMS Open Data with a trigger 

JPSIAnaOpenData2011_Model_Pythia.py
This file work with "JpsiAnalyzerOpen2011_2.cc". This script in python was made to get acess the online root files in CMS Open Data withot a trigger 

loop_trigger_normal.py
This script in PYTHON was made to work with the files above. With this you can chose the type of data you will acess anda have outputs for each online file.root. That way you can automate the acess in each adress, besides you can know what file/adress has stop if stops.

analysis.cc
Program to extract data (variables and vector format) of the root file and create histograms. This root file is made of others root files, therefore, has multiple entries which we access with a diferent method than a normal root file.
