import os
import subprocess
import pandas as pd
from io import StringIO
import numpy as np
import argparse
import sys

"""Este script es para quedarnos con los picos solapantes entre control y tratado ( u otro par de archivos). Con una base que solape es suficiente. Luego tambien permite calulcar el socre y FE mediante el promedio de los score de los picos que definen un pico consenso. Finalmente calculo el numero de picos comunes entre los samples y los picos unicos"""


def create_command(gsize,peakscondition1,peakscondition2):
    samples =  peakscondition1+ " " + peakscondition2
    command = "/software/HOMER/bin/mergePeaks -d given -gsize " + str(gsize) + " " + samples
    return command

def runcml(cml): ##runnea el comando creado en la anterior definicion
    """Runs a command"""
    process = subprocess.Popen(cml.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout, stderr = process.communicate()
    stderr = stderr.decode("utf-8")
    stdout = stdout.decode("utf-8")
    df = pd.read_table(StringIO(stdout), sep="\t",skiprows=1,header=None,names=["PeakName","chromosome","start","end","strand","score","samples","npeaks","sample1","sample2"])
    df = df.sort_values(["chromosome","start","end"],ascending=[True,True,True])
    print(df.head())
    return df

def create_sample_df(peakscondition1,peakscondition2): #crear la nueva base de datos con los dos tratamientos
    colnames = ["chromosome","start","end","PeakName","score","strand","FE","qValue","replicates","samples_ID","npeaks"]
    df = pd.DataFrame([], columns=colnames)
    df = pd.concat([df,pd.read_table(peakscondition1, sep="\t",names=colnames),pd.read_table(peakscondition2, sep="\t",names=colnames)])
    #juntar en una base las dos bases de datos que teniamos
    df = df.sort_values(["chromosome","start","end"],ascending=[True,True,True])
    return df

def create_mergedPeaks(df,df_indv):
    rows = []
    countPeak = 0
    df = df.fillna("-") #sustituir los valores perdidos con '-'
    for index,row in df.iterrows():
        countPeak += 1
        PeakName = "Peak" + str(countPeak)
        replicates = len(row["samples"].split("|"))
        values_ = row[["sample1","sample2"]].values.tolist()
        values_ = [x for x in values_ if x != "-" and str(x) != 'nan']
        values_ = ",".join(values_).split(",")
        samples_ID = [x.split("_")[1] for x in values_] 
        samples_ID.sort()
        samples_ID = ",".join(set(samples_ID))
        score = int(np.median(df_indv[df_indv["PeakName"].isin(values_)]["score"].tolist()))
        try:
            score_control = int(np.median(df_indv[(df_indv["PeakName"].isin(values_)) &(df_indv["Treatment"]=="C")]["score"].tolist()))
        except:
            score_control = 0
        try:
            score_treatment = int(np.median(df_indv[(df_indv["PeakName"].isin(values_)) &(df_indv["Treatment"]=="P")]["score"].tolist()))
        except:
            score_treatment = 0
        try:
            FE_control = (np.median(df_indv[(df_indv["PeakName"].isin(values_)) &(df_indv["Treatment"]=="C")]["FE"].tolist()))
            FE_control = round(FE_control,2)
        except:
            FE_control = 0
        try:
            FE_treatment = (np.median(df_indv[(df_indv["PeakName"].isin(values_)) &(df_indv["Treatment"]=="P")]["FE"].tolist()))
            FE_treatment = round(FE_treatment,2)
        except:
            FE_treatment = 0
        if score_control !=0 and score_treatment != 0:#si ninguno es 0, calculamos el ratio
            scorerate = float(round(np.log(score_treatment/score_control),3)) # treatment vs control
        else:
            scorerate = 0
        if len(df_indv[(df_indv["PeakName"].isin(values_)) & (df_indv["Treatment"]=="C")]) == 0:
            replicates_Control = 0
        else:
            replicates_Control = len(set(",".join(df_indv[(df_indv["PeakName"].isin(values_)) & (df_indv["Treatment"]=="C")]["samples_ID"].tolist()).split(",")))
        if len(df_indv[(df_indv["PeakName"].isin(values_)) & (df_indv["Treatment"]=="P")]) == 0:
            replicates_Treated = 0
        else:
            replicates_Treated = len(set(",".join(df_indv[(df_indv["PeakName"].isin(values_)) & (df_indv["Treatment"]=="P")]["samples_ID"].tolist()).split(",")))
        chromosome = row["chromosome"]
        start = row["start"]
        end = row["end"]
        npeaks = row["npeaks"]
        strand = "."
        if row["sample1"] != "-":
            sample1 = "YES"
        else:
            sample1 = "NO"
        if row["sample2"] != "-":
            sample2 = "YES"
        else:
            sample2 = "NO"
        rows.append([chromosome,start,end,PeakName,score,score_control,score_treatment,scorerate,FE_control,FE_treatment,strand,replicates,samples_ID,npeaks,replicates_Control,replicates_Treated,sample1,sample2])
    return pd.DataFrame(rows,columns=["chromosome","start","end","PeakName","score","scoreControl","scoreTreatment","scoreRate","FE_control","FE_treatment","strand","replicates","samples_ID","npeaks","replicates_Control","replicates_Treated","sample1","sample2"])

def define_areas(df_mergedPeaks):
    area1 = len(df_mergedPeaks[df_mergedPeaks["sample1"] == "YES"]) # All sample1
    area2 = len(df_mergedPeaks[df_mergedPeaks["sample2"] == "YES"]) # All sample2
    n12 = len(df_mergedPeaks[(df_mergedPeaks["sample1"] == "YES") & (df_mergedPeaks["sample2"] == "YES")]) #  sample1 n sample2 (intersection)
    return area1,area2,n12


gsize = "1.87e9"
#sex = sys.argv[1] #sys.argv es una lista en Python, que contiene los argumentos de la línea de comandos pasados al script.
sextype = "M" #"F"
peakscondition1 = "C" + sextype+ "_HomerMergedPeaks.bed"
peakscondition2 = "P"+ sextype+ "_HomerMergedPeaks.bed"
samples = [peakscondition1.split("_")[0], peakscondition2.split("_")[0]]
outputdir = "/media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/peaks/blacklist/"

#Dos tipos full o partial (full es en todas las replicas y partial en 2 o mas)

cml = create_command(gsize,peakscondition1,peakscondition2)
df = runcml(cml)
df_indv = create_sample_df(peakscondition1,peakscondition2)
df_indv["Treatment"] = [x.split("_")[0][0] for x in df_indv.PeakName.tolist()]
df_mergedPeaks = create_mergedPeaks(df,df_indv)


df_mergedPeaks["FERate"] = df_mergedPeaks["FE_treatment"]/df_mergedPeaks["FE_control"]
df_mergedPeaks["FERate"] = round(df_mergedPeaks["FERate"],2)


with open(sextype + "_venndiagram.log","w") as fo:
     area1,area2,n12 = define_areas(df_mergedPeaks)
     header = "\t".join(map(str,["area1","area2","n12"])) + "\n"
     line = "\t".join(map(str,[area1,area2,n12]))
     fo.write(header)
     fo.write(line)


try:
    os.makedirs("IGV") #El método os.makedirs() en Python se usa para crear un directorio de forma recursiva
except:
    pass

df_mergedPeaks = df_mergedPeaks[['chromosome', 'start', 'end', 'PeakName', 'score', 'scoreControl','scoreTreatment', 'scoreRate', 'FE_control', 'FE_treatment',"FERate",'samples_ID','replicates_Control','replicates_Treated']]
df_mergedPeaks = df_mergedPeaks.fillna(0) ##rellenar los valores perdidos con un 0
##### Filter

df = df_mergedPeaks

common = df[(df["replicates_Control"]>= 2) & (df["replicates_Treated"]>=2)] #el archivo que contiene los picos de control
treated = df[(df["replicates_Control"]== 0) & (df["replicates_Treated"]>=2)]
control = df[(df["replicates_Control"]>= 2) & (df["replicates_Treated"]==0)]


common["Type"] = "Common"
treated["Type"] = "Treated"
control["Type"] = "Control"


df = pd.concat([common,control,treated])

df.to_csv( "Peaks.bed",sep="\t",header=False,index=False)
common.to_csv( "IGV/Common.bed",sep="\t",header=False,index=False)
control.to_csv("IGV/Control.bed",sep="\t",header=False,index=False)
treated.to_csv("IGV/Treated.bed",sep="\t",header=False,index=False)
