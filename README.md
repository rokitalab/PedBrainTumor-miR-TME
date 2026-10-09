# Pan-pediatric brain tumor microRNA networks define immune microenvironment states and responses to CAR T-cell therapy

Kaleem Coleman^, Ryan Corbett^, Zhongzhen Yi, Bicna Song, Alex Sickler, Aylar Babaei, Jo Lynne Rokita+, Dalia Haydar+*

^ equal authorship
+ co-senior authors
* correspondence: Dalia Haydar, dhaydar@childrensnational.org

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses

1. Clone the repository
```
git clone git@github.com:rokitalab/PedBrainTumor-miR-TME.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/rokita-lab/haydar-mirna:v1.0.3
```
NOTE: if running on a Mac with Apple Silicon chip (M1-M4), please add `--platform linux/amd64`; otherwise add `--platform linux/arm64`

3. Start the docker container, from the `PedBrainTumor-miR-TME` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/haydar-mirna pgc-images.sbgenomics.com/rokita-lab/haydar-mirna:v1.0.3
```
NOTE: if running on a Mac with Apple Silicon chip (M1-M4), please add `platform linux/amd64`

4. To execute shell within the docker image, from the `PedBrainTumor-miR-TME` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the `download-data.sh` shell script to obtain latest data files: 
```
bash download_data.sh
```

6. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/PedBrainTumor-miR-TME/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── LICENSE
├── README.md
├── analyses
├── data
├── download_data.sh
├── figures
└── scripts
```

## Miniconda End User License Agreement

Copyright Notice: Miniconda® © 2015, Anaconda, Inc. 

All rights reserved. Miniconda® is licensed, not sold. 

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer;

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution; 

3. The name Anaconda, Inc. or Miniconda® may not be used to endorse or promote products derived from this software without specific prior written permission from Anaconda, Inc.; and 

4. Miniconda® may not be used to access or allow third parties to access Anaconda package repositories if such use would circumvent paid licensing requirements or is otherwise restricted by the Anaconda Terms of Service.

DISCLAIMER: THIS SOFTWARE IS PROVIDED BY ANACONDA “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE , AND NON-INFRINGEMENT ARE DISCLAIMED. IN NO EVENT SHALL ANACONDA BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF MINICONDA®, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Code Authors

Ryan Corbett ([@rjcorb](https://github.com/rjcorb)), Bicna Song ([@bicnasong](https://github.com/bicnasong))
