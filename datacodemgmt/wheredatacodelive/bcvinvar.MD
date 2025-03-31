---
geometry: margin = 30mm
fontsize: 12pt
---

__Date of Study__  
May 2024 - August 2024

__Date of Public Archiving:__  
End of 2024 (hopefully, when wine-growing area shapefiles are published by Nacho and his student)

__Last Modified:__  
31 March 2025

__Goal__  

Assess suitability of common winegrape varieties in wine-growing areas of BC using models based on climatic limits defined by phenological and climate data for 1994 - 2023

__Contributors__

Ken Michiko Samson

Elizabeth Wolkovich

__General Files__

| __File__  |__Where__ |__What__  |
|:---       |:-----    |:----     |
| BCWGC Report (full and short version) | [Full](https://stateofwine.org/wp-content/uploads/2024/09/bcvarsuitability.pdf) and [snapshot](https://stateofwine.org/wp-content/uploads/2024/09/bcvarsuitability_short.pdf) reports available online at stateofwine.org; main working files in the private GitHub repo [bcvinvarieties](https://github.com/lizzieinvancouver/bcvinvarieties): docs/ | Report for the British Columbia Wine Grape Council |

__Data and Code__

Main workspace in GitHub repo bcvinvar: https://github.com/lizzieinvancouver/bcvinvarieties
Big data files (indicated separately below) in midge

| __File__  |__Where__ |__What__  |
|:---       |:-----    |:----     |
| Data | GitHub: analyses/input | Phenological data from vin, wine-growing region shapefiles from Brad and Nacho, and climate data from EnvCanada|
| Gridded Climate Data | midge: analyses/input/ | Historical climate data from Daymet and projected climate data from NEX-GDDP-CMIP6 dataset |
| Data Formatting Code | GitHub: analyses/data_formatting/ | Code to format data |
| Data Preparation Code | GitHub: analyses/data_preparation/ | Code to run models |
| Formatted Data | GitHub: analyses/output/formatted_data/ | Formatted data used by data preparation scripts |
| Model Data | GitHub: analyses/output/prepared_data/ | Model data and projections |
| General Climate Envelope Model Data | midge: analyses/output/prepared_data/ | Model output and metrics of general climate envelope model  |
| Winegrape Climate Envelope Data | midge: analyses/output/prepared_data/wr_data/ | Model output of winegrape climate envelope model |
| Figure-making code | GitHub: analyses/ | Code to make report figures |
| Figures | GitHub: analyses/figures | Figures for report |
| Hardiness Models | GitHub: analyses/hardiness/ | Hardiness models used by Faith and Al |
| Checks | GitHub: analyses/checks/ | Code used to check data and test model-running scripts |

__Notes__

Additional contributions are made by:

- Faith Jones and Al Kovaleski shared hardiness models and data and assisted with cold hardiness modeling

- Ignacio Morales-Castilla shared shapefiles of wine-growing regions and assisted with species distribution modeling

- Justin Ngo made Figure 1 for the report

- Britany Wu reviewed literature on winegrowing in British Columbia

- Carl Bogdanoff

- Ben-Min Chang shared shapefiles of wine-growing regions courtesy of Scott Smith

- Brad Estergaard  shared shapefiles of wine-growing regions

- Steve Marsh