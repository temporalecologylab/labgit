_Date of study_  
Start date (**?**) - today

_Date of Public Archiving:_  
Enter when data should be (or was) publicly-archived following lab data management plan.  (two years after completion of project)

_Last Modified:_  
25 February 2025

_Goal_  
Quantitative review of literature on the response of seed germination to environmental cues.

_Contributors_  
**To be completed:**  
Lizzie Wolkovich - [e.wolkovich@ubc.ca](mailto:e.wolkovich@ubc.ca)    
Daniel Buonaiuto - [dbuonaiuto@umass.edu](mailto:dbuonaiuto@umass.edu)  
Deirdre Loughnan - [deirdre.loughnan@gmail.com](mailto:deirdre.loughnan@gmail.com)  
Xiaomao Wang - [wangxm29@mail.ubc.ca](mailto:wangxm29@mail.ubc.ca)  
Christophe Rouleau-Desrochers - [christophe.rouleaudesrochers@ubc.ca](mailto:christophe.rouleaudesrochers@ubc.ca)  
Justin Ngo - [ngojustintwz1@gmail.com](mailto:ngojustintwz1@gmail.com)  
Ken Michiko Samson - [shmkkmhs12@gmail.com](mailto:shmkkmhs12@gmail.com)  
Frederik Baugmarten - [frederik.baumgarten@wsl.ch ](mailto:frederik.baumgarten@wsl.ch)  
Bunch of people...  
Victor Van der Meersch - [victor.vandermeersch@cefe.cnrs.fr](mailto:victor.vandermeersch@cefe.cnrs.fr)  

_General Files_

| __File__  |__Where__ |__What__  | 
|:---       |:---      |:---      |  
| Overview | [GitHub: readme](https://github.com/lizzieinvancouver/egret/blob/main/_README.md) | Overview of general goals and update |
| General methods | [GitHub: notes](https://github.com/lizzieinvancouver/egret/tree/main/notes) | Info on how to start a meta-analysis and scrape data |

_Data and Code_

__Github__  [https://github.com/lizzieinvancouver/egret](https://github.com/lizzieinvancouver/egret)  
__W: drive__ Weldshare > Wolkovich Lab > ?

| __File__  |__Where__ |__What__  |
|:---       |:---      |:---      |
| Raw dataset `egretData.csv` | [GitHub: analyses/input](https://github.com/lizzieinvancouver/egret/blob/main/analyses/input/egretData.csv) | Raw scrapped data from published experiments |
| Cleaning files | [GitHub: analyses/cleaning](https://github.com/lizzieinvancouver/egret/blob/main/analyses/input/egretData.csv) | Set of scripts to clean data (from published experiments), gathered in `cleanall.R` (pipeline to define how and in what order they should be run) |
| Raw USDA dataset | [GitHub: analyses/scrapeUSDAseedmanual/input](https://github.com/lizzieinvancouver/egret/blob/main/analyses/scrapeUSDAseedmanual/input) | Raw scrapped tables from USDA seed manual |
| Cleaning files for USDA data | [GitHub: analyses/scrapeUSDAseedmanual/cleaning](https://github.com/lizzieinvancouver/egret/blob/main/analyses/scrapeUSDAseedmanual/cleaning) | Set of scripts to clean USDA data, gathered in `cleanAllUSDA.R` |
| Cleaning misc. files | [GitHub: analyses/cleaningOther](https://github.com/lizzieinvancouver/egret/blob/main/analyses/cleaningOther) | Set of scripts to clean other data, including Baskin database |
| Data filtering and selection | [GitHub: analyses/studyDesign](https://github.com/lizzieinvancouver/egret/blob/main/analyses/studyDesign) | Set of scripts to compare treatments across studies and develop a framework for selecting the most relevant combinations (`decisionRules.R`) |
| Exploratory work on models | [GitHub: analyses/modeling](https://github.com/lizzieinvancouver/egret/blob/main/analyses/modeling) | Set of scripts to compare Binomial/Beta distribution with simulated data, and fit the very first models to real data |
| Stan models | [GitHub: analyses/stan](https://github.com/lizzieinvancouver/egret/blob/main/analyses/stan) | Stan models |

__Possible Extras:__

Any amendments to when public archiving happen should mentioned here and an asterisk given above where archiving date is given.  

Be sure all your data is somewhere where it is backed up as per the data management plan.  

Check this file for accuracy, and update as needed, every 6 months or sooner.  
  
  
  
  
  
![CC License](ccimage/CC_BY-SA_icon.svg.png)
