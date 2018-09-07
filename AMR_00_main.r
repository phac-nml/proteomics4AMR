#Copyright Government of Canada 2018
#
#"Written by: Julie Chih-yu Chen @ National Microbiology Laboratory, Public Health Agency of Canada"
#
#"Licensed under the Apache License, Version 2.0 (the ""License""); you may not use"
#this work except in compliance with the License. You may obtain a copy of the
#License at:
#
#http://www.apache.org/licenses/LICENSE-2.0
#
#"Unless required by applicable law or agreed to in writing, software distributed"
#"under the License is distributed on an ""AS IS"" BASIS, WITHOUT WARRANTIES OR"
#"CONDITIONS OF ANY KIND, either express or implied. See the License for the"
#specific language governing permissions and limitations under the License.


setwd("D:\\ClientData\\_AMR\\_scriptFinal") 

### summarize RGI results
source("AMR_01_RGI.r")

### summarize ResFinder results
source("AMR_02_ResFinder.r")

### Analyze proteomics results from MaxQuant
source("AMR_03_functions.r")
source("AMR_04_proteomics.r")


