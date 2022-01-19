import pandas as pd
from collections import defaultdict
from itertools import combinations
from os.path import join
import os


#====config and sample=====
configfile: "config.yaml"

def getList(partnumber):
    LIST = []
    for i in range(1,partnumber+1):
            j = str(i)
            LIST.append(j.rjust(3,'0'))
    return LIST

outdir = config["outdir"]
Isamples = config["rawdata"]["ISO"]
Rsamples = config["rawdata"]["RNA"]
# PEsample = Ssamples + Rsamples + Hsamples
# PEdir=[]
# for i in Ssamples :
#     PEdir[i]="01.survey"
# for i in Rsamples :
#     PEdir[i]="01.survey"

name = config["name"]

species = config["parameters"]["EDTA"]["species"]
step = config["parameters"]["EDTA"]["step"]
sensitive = config["parameters"]["EDTA"]["sensitive"]
lib = config["parameters"]["EDTA"]["lib"]
cds = config["parameters"]["EDTA"]["cds"]

# samples = pd.read_csv(config["samples"],sep="\t").set_index("sample",drop=False)
# CHR = config["Chr"]
# TYPE = config["type"]
# Title = config["title"]
# POP = config["pop"]
#contigs = pd.read_csv(config["ref"]["genome"] + ".fai",header=None, usecols=[0], squeeze=True, dtype=str,sep="\t")




#========biosoft========

module = config["biosoft"]["env"]["module"]
mf = config["biosoft"]["env"]["mf"]
mp = config["biosoft"]["env"]["mp"]
asm = config["biosoft"]["env"]["asm"]
#env
EDTA = config["biosoft"]["repeat"]["EDTA"]


# braker
mapping = config["biosoft"]["braker"]["mapping"]
QC = config["biosoft"]["braker"]["QC"]
braker2 = config["biosoft"]["braker"]["braker2"]
genemark = config["biosoft"]["braker"]["genemark"]
augustus_conf = os.getcwd()+"/"+config["biosoft"]["braker"]["augustus_conf"]



maker = config["biosoft"]["maker"]["maker2"]

