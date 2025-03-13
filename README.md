# ChemEng277-Project
A PK Data predictor for small molecule oral dose

Hello and welcome to my pharmacokinetic (PK) data model! The objective of this code is to predict the PK Data of a
small molecule orally delivered in a human. This is done using specific properties of the drug established by the user
and clustering if inputted by the user. Current data inputted is the data used for the accompanying paper and the Trials
folder holds a collection of data used for analysis in the said paper. Please refer to the "REQUIRMENTS.txt" file for
required modules for this model.

To use the model it is important to go over the specific items used within this package:
-"TrainData" Folder:
    This is where the time course PK graphs downloaded as json files from PK-DB.com are inputted into the model.
    The files must be named using this convention: DRUGNAME_DOSE.json, where dose is in mg.
    Additionally, the "Drug_Properties.xlsx" file is located here. This is an excel file where the desired properties to
    be compared are stored in the following column format: DRUGNAME ... Property 1 ... Property n. DRUGNAMES must be
    organized alphabetically.
-"TestData" Folder:
    This is where the time course PK graphs downloaded as json files from PK-DB.com are inputted into the model.
    The files must be named using this convention: DRUGNAME_DOSE.json, where dose is in mg. The files must be named
    using this convention: DRUGNAME_DOSE.json, where dose is in mg.
    Additionally, the "Test_Drug_Properties.xlsx" file is located here. This is an excel file where the desired
    properties to be compared are stored in the following column format: DRUGNAME ... Property 1 ... Property n ...
    PREDICTED DOSE. DRUGNAMES must be organized alphabetically and only one dose predicted at a time.
-"Results_TrainData" Folder:
    This folder is populated with excel files titled "DRUGNAME.xlsx" containing the actual and predicted PK values
    for that drug used to generate the model. The cluster groups and mean property values are listed
    in "cluster_properties.xlsx".
-"Results_TestData" Folder:
    This folder is populated with excel files titled "DRUGNAME.xlsx" containing the actual and predicted PK values
    for that drug predicted from the model.

The program takes in 6 arguments:
    1. file path to the training data folder ("TrainData")
    2. file path to the test data folder ("TestData")
    3. file path to the training data results folder ("Results_TrainData")
    4. file path to the test data results folder ("Results_TestData")
    5. int value for time cutoff on PK data (can be set to 1000 to avoid cutoff)
    6. int value for number of drug clusters used in modeling

example of running the model for a 90 minutes cutoff and 5 clusters:
.\CHEME277Project.py
C:\Users\goehr\OneDrive\Desktop\Stanford\ChemE277\FinalProject\TrainData
C:\Users\goehr\OneDrive\Desktop\Stanford\ChemE277\FinalProject\TestData
C:\Users\goehr\OneDrive\Desktop\Stanford\ChemE277\FinalProject\Results_TrainData
C:\Users\goehr\OneDrive\Desktop\Stanford\ChemE277\FinalProject\Results_TestData
90
5
