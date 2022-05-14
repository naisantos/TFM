import os
import subprocess
import pandas as pd
from io import StringIO
import numpy as np
import argparse

"""Este script es para quedarnos con los picos solapantes entre replicados biologicos. Con una base que solape es suficiente. Luego tambien permite calulcar el socre mediante el promedio de los score de los picos que definen un pico consenso. Finalmente calculo el numero de picos comunes entre los samples y los picos unicos"""

def create_command(call_peak_dir,sampletype,gsize):
    samples = [x for x in os.listdir(call_peak_dir) if sampletype in x and ".noBlackLists.narrowPeak" in x]
    nsample = len(samples)
    samples.sort()
    samples = " ".join(samples)
    command = "/software/HOMER/bin/mergePeaks -d given -gsize " + str(gsize) + " " + samples 
    #los arguementos necesarios que se utilizan para HOMER
    return command,nsample

def runcml(cml,nsample):
    """Runs a command"""
    process = subprocess.Popen(cml.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout, stderr = process.communicate()
    stderr = stderr.decode("utf-8")
    stdout = stdout.decode("utf-8")
    if nsample ==2:
        df = pd.read_table(StringIO(stdout), sep="\t",skiprows=1,header=None,names=["PeakName","chromosome","start","end","strand","score","samples","npeaks","sample1","sample2"])
    else:
        df = pd.read_table(StringIO(stdout), sep="\t",skiprows=1,header=None,names=["PeakName","chromosome","start","end","strand","score","samples","npeaks","sample1","sample2","sample3"])
    df = df.sort_values(["chromosome","start","end"],ascending=[True,True,True])
    return df

# HOMER FILE
# 1. Merged Peak name (will start with "Merged-")
# 2. chromosome
# 3. start (average from merged peaks)
# 4. end (average from merged peaks)
# 5. strand
# 6. Average peak score (actually, the average of the original values in column 6 of the peak files - or column 5 of BED files)
# 7. Original peak files contributing to the merged peak
# 8. Total number of peaks merged (occasionally more than one peak from a single file will be merged if the peaks are within the specify distance or two or more peaks from one file overlap with the same single peak(s) from another file)

def create_sample_df(call_peak_dir,sampletype): #creamos un solo archivo para las 3 replicas
    samples = [x for x in os.listdir(call_peak_dir) if sampletype in x and ".noBlackLists.narrowPeak" in x]
    nsample = len(samples)
    samples.sort()
    colnames = ["chromosome","start","end","name","score","strand","signalValue","pvalue","qvalue","peak"]
    df = pd.DataFrame([], columns=colnames)
    if nsample == 2:
        df = pd.concat([df,pd.read_table(samples[0], sep="\t",names=colnames),pd.read_table(samples[1], sep="\t",names=colnames)])
        df["pvalue"] = 10**-df["pvalue"]
        df["qvalue"] = 10**-df["qvalue"]
    elif nsample == 3:
        df = pd.concat([df,pd.read_table(samples[0], sep="\t",names=colnames),pd.read_table(samples[1], sep="\t",names=colnames),pd.read_table(samples[2], sep="\t",names=colnames)])
        df["pvalue"] = 10**-df["pvalue"]
        df["qvalue"] = 10**-df["qvalue"]
    df = df.sort_values(["chromosome","start","end"],ascending=[True,True,True])
    return df


def get_entichment(peakidlist):
    fe,qvalue =df_indv[df_indv["name"].isin(peakidlist)].sort_values("score",ascending=False)[["signalValue","qvalue"]].iloc[0].tolist()
    return fe,qvalue

def create_mergedPeaks(df,df_indv,nsample,sampletype):
    rows = []
    countPeak = 0
    df = df.fillna("-")
    for index,row in df.iterrows(): #recorremos cada archivo 
        countPeak += 1
        PeakName = sampletype + "_Peak" + str(countPeak)
        replicates = len(row["samples"].split("|"))
        if nsample == 3:
            values_ = row[["sample1","sample2","sample3"]].values.tolist()
        elif nsample == 2:
            values_ = row[["sample1","sample2"]].values.tolist()
        values_ = [x for x in values_ if x != "-"]
        values_ = ",".join(values_).split(",")
        samples_ID = [x.split("_")[1] for x in values_]
        samples_ID.sort()
        samples_ID = ",".join(set(samples_ID))
        score = int(np.mean(df_indv[df_indv["name"].isin(values_)]["score"].tolist()))
        fe,qvalue = get_entichment(values_)
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
        if nsample == 3:
            if row["sample3"] != "-":
                sample3 = "YES"
            else:
                sample3 = "NO"
            rows.append([chromosome,start,end,PeakName,score,strand,fe,qvalue,replicates,samples_ID,npeaks,sample1,sample2,sample3])
        else:
            rows.append([chromosome,start,end,PeakName,score,strand,fe,qvalue,replicates,samples_ID,npeaks,sample1,sample2])
    if nsample == 2:
        return pd.DataFrame(rows,columns=["chromosome","start","end","PeakName","score","strand","FE","qValue","replicates","samples_ID","npeaks","sample1","sample2"])
    elif nsample == 3:
        return pd.DataFrame(rows,columns=["chromosome","start","end","PeakName","score","strand","FE","qValue","replicates","samples_ID","npeaks","sample1","sample2","sample3"])


parser = argparse.ArgumentParser()
parser.add_argument("sampletype", type=str,help="prefix to call the sample type") #en nuestro caso, pueden ser CM, CF, PM, PF
parser.add_argument("gsize", type=str, default="1.87e9", help="genome size. Can be exponential notation. Default mm39 size")
parser.add_argument("call_peak_dir", type=str, default=".",help="Directory where the blacklisted narrow peaks are located")
parser.add_argument("-o", "--output", action="store", default="HomerMergedPeaks.bed")
args = parser.parse_args()


cml,nsample = create_command(args.call_peak_dir,args.sampletype,args.gsize)

df_indv = create_sample_df(args.call_peak_dir,args.sampletype)

df = runcml(cml,nsample) # The score calulcated by HOMER is not functioning so we add an extra function for that

df_mergedPeaks = create_mergedPeaks(df,df_indv,nsample,args.sampletype)




def define_areas(nsample,df_mergedPeaks):
    area1 = len(df_mergedPeaks[df_mergedPeaks["sample1"] == "YES"]) # All sample1
    area2 = len(df_mergedPeaks[df_mergedPeaks["sample2"] == "YES"]) # All sample2
    if nsample == 3:
        area3 = len(df_mergedPeaks[df_mergedPeaks["sample3"] == "YES"]) # All sample3
        n12 = len(df_mergedPeaks[(df_mergedPeaks["sample1"] == "YES") & (df_mergedPeaks["sample2"] == "YES")])  #  sample1 n sample2 (intersection)
        n13 = len(df_mergedPeaks[(df_mergedPeaks["sample1"] == "YES")  & (df_mergedPeaks["sample3"] == "YES")])  #  sample1 n sample3 (intersection)
        n23 = len(df_mergedPeaks[(df_mergedPeaks["sample2"] == "YES") & (df_mergedPeaks["sample3"] == "YES")])  #  sample2 n sample3 (intersection)
        n123 = len(df_mergedPeaks[(df_mergedPeaks["sample1"] == "YES") & (df_mergedPeaks["sample2"] == "YES") & (df_mergedPeaks["sample3"] == "YES")])  #  sample1 n sample2 n sample3 (intersection)
        return area1,area2,area3,n12,n13,n23,n123
    else:
        n12 = len(df_mergedPeaks[(df_mergedPeaks["sample1"] == "YES") & (df_mergedPeaks["sample2"] == "YES")]) #  sample1 n sample2 (intersection)
        return area1,area2,n12

with open("venndiagram_" + args.sampletype + ".log","w") as fo: #creamos el archivo donde aparecerá el número de picos de cada replica,
    #y cuántos se solapan entre ellos
    if nsample == 3:
        area1,area2,area3,n12,n13,n23,n123 = define_areas(nsample,df_mergedPeaks)
        header = "\t".join(map(str,["area1","area2","area3","n12","n13","n23","n123"])) + "\n"
        line = "\t".join(map(str,[area1,area2,area3,n12,n13,n23,n123]))
    elif nsample == 2:
        area1,area2,n12 = define_areas(nsample,df_mergedPeaks)
        header = "\t".join(map(str,["area1","area2","n12"])) + "\n"
        line = "\t".join(map(str,[area1,area2,n12]))
    fo.write(header)
    fo.write(line)

if nsample==3:
        df_mergedPeaks[df_mergedPeaks.columns.tolist()[:-3]].to_csv(args.output,sep="\t",header=False,index=False) #creamos el nuevo archivo

else:
        df_mergedPeaks[df_mergedPeaks.columns.tolist()[:-2]].to_csv(args.output,sep="\t",header=False,index=False)


# python3 HomerMergedPeaks.py -o CM_HomerMergedPeaks.bed CM 1.87e9 /media/sequentia/synology_office/Naiara/H3K4me3/bam_files/Filtered/peaks/blacklist
