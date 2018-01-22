#!/usr/bin/env python3
import mygene

import numpy as np
import scipy as sp
import pandas as pd
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import tempfile
import os
from rpy2.robjects import pandas2ri
pandas2ri.activate()
import string
import rpy2
from rpy2 import robjects


functions = '''read.peptides <- function(dat, cha){
  output <- NULL

  dat$Sequence <- as.character(dat$Sequence)
  dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
  dat$Quan.Usage <- as.character(dat$Quan.Usage)
  dat$Quan.Info <- as.character(dat$Quan.Info)
  dat$Isolation.Interference <- as.numeric(as.character(dat$Isolation.Interference))

  dat <- subset(dat, Isolation.Interference<=30)
  dat <- subset(dat, Quan.Usage=="Used")
  dat <- subset(dat, Protein.Group.Accessions!="")
  dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))
}

quantify.proteins <- function(dat, cha){
  e.function <- function(x, seq) tapply(x, seq, median)
  output <- NULL

  dat$Sequence <- toupper(dat$Sequence) # Capital letters
  accessions <- as.character(unique(dat$Protein.Group.Accessions))
  n.proteins <- length(accessions)
  n.cha <- length(cha)

  for(k in 1:n.proteins){
    id <- accessions[k]
    sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
    sdat[cha] <- log2(sdat[cha])
    sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
    pdat <- sdat[, -1]
    n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
    temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])
    n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)
    if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
    pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
    output <- rbind(output, pdat)
  }
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- round(output[,1:n.cha],3)
  row.names(output) <- accessions
  output <- as.data.frame(output)

  return(output)
}

eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

trt.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- topTreat(treat(fit), coef=2,lfc=1)
  return(fit.eb)
}

eb.fit.mult <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coef[, "tr2"]
  df.0 <- rep(fit.eb$df.prior, n)
  df.r <- fit.eb$df.residual
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
  t.mod <- fit.eb$t[, "tr2"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "tr2"]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb.mult <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.0, df.r, s2.0, s2, s2.post)
  results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
  return(results.eb.mult)
}'''

def boxplot(df, columns , xval, yval, png): # groups is a dict to group values by
    df_ = df[columns]
    rdf =pandas2ri.py2ri(df_)
    ro.globalenv['data'] = rdf
    c = "png('{}',width = 12, height = 8, units = 'in', res = 300)".format(png); ro.r(c)
    c = 'print(boxplot(data, las = 2)); dev.off()'; ro.r(c)

def list_density(data, title, outpath):
    rdf =pandas2ri.py2ri(data)
    ro.globalenv['data'] = rdf
    c = "attach(data)"; ro.r(c)
    c = "library(ggplot2)"; ro.r(c)
    c = "p <- ggplot(aes(x=PEP, colour=Type), data=data)"; ro.r(c)
    c = "pl1 <- p + geom_density()+xlab('Posterior error probability (PEP) score')+ggtitle('{}')".format(title); ro.r(c)
    c = "png('{}',width = 12, height = 8, units = 'in', res = 300)".format(outpath); ro.r(c)
    c = "print(plot(pl1)); dev.off()"; ro.r(c)

def qval(dataframe, pval):   # pandas df, treated columns, control columns
    df = dataframe.copy()
    #col_dct = []
    #count = 0
    #for col in df:
    #    col_dct[count] = col
    #    df.rename(columns={col : count}, inplace =True)
    #    count += 1
    c="library(qvalue)"; ro.r(c)
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c="data$q.value <- qvalue(data${})$q".format(pval); ro.r(c)
    c = 'data'
    dfebi = pandas2ri.ri2py(ro.r[c])
    #dfebi.rename(columns = col_dct, inplace = True)
    return dfebi

def spearman_plot(df_, col1, col2, path):   # pandas df, treated columns, control columns
    df = df_.copy()
    exclude = set(string.punctuation)
    col1_ = ''.join(ch for ch in col1 if ch not in exclude)
    col2_ = ''.join(ch for ch in col2 if ch not in exclude)
    df.rename(columns={col1:col1_,col2:col2_},inplace=True)

    c="library(qvalue)"; ro.r(c)
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c = "attach(data)"; ro.r(c)
    c = "my.cor=cor.test({},{},method='spearman',use='pairwise.complete.obs')".format(col1_, col2_);ro.r(c)
    c = "png('{}')".format(path+'spearman.png'); ro.r(c)
    c = "print(plot({}, {}, xlab='{}', ylab='{}', pch=21)); dev.off()".format(col1_, col2_, col1_,col2_,col2_,col1_); ro.r(c)


def spearman(df_, col1, col2):   # pandas df, treated columns, control columns
    df = df_.copy()
    exclude = set(string.punctuation)
    col1_ = ''.join(ch for ch in col1 if ch not in exclude)
    col2_ = ''.join(ch for ch in col2 if ch not in exclude)
    df.rename(columns={col1:col1_,col2:col2_},inplace=True)

    c="library(qvalue)"; ro.r(c)
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c = "attach(data)"; ro.r(c)
    c = "my.cor=cor.test({},{},method='spearman',use='pairwise.complete.obs')".format(col1_, col2_);ro.r(c)
    pval = float(rpy2.robjects.r('my.cor$p.value')[0])
    estimate = float(rpy2.robjects.r('my.cor$estimate')[0])
    return estimate, pval

def pearson_plot(df_, col1, col2, path):   # pandas df, treated columns, control columns
    df = df_.copy()
    exclude = set(string.punctuation)
    col1_ = ''.join(ch for ch in col1 if ch not in exclude)
    col2_ = ''.join(ch for ch in col2 if ch not in exclude)
    df.rename(columns={col1:col1_,col2:col2_},inplace=True)

    c="library(qvalue)"; ro.r(c)
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c = "attach(data)"; ro.r(c)
    c = "my.cor=cor.test({},{},method='pearson',use='pairwise.complete.obs')".format(col1_, col2_);ro.r(c)
    c = "png('{}')".format(path+'pearson.png'); ro.r(c)
    c = "print(plot({}, {}, xlab='{}', ylab='{}', pch=21)); dev.off()".format(col1_, col2_, col1_,col2_,col2_,col1_); ro.r(c)


def pearson(df_, col1, col2):   # pandas df, treated columns, control columns
    df = df_.copy()
    exclude = set(string.punctuation)
    col1_ = ''.join(ch for ch in col1 if ch not in exclude)
    col2_ = ''.join(ch for ch in col2 if ch not in exclude)
    df.rename(columns={col1:col1_,col2:col2_},inplace=True)

    c="library(qvalue)"; ro.r(c)
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c = "attach(data)"; ro.r(c)
    c = "my.cor=cor.test({},{},method='pearson',use='pairwise.complete.obs')".format(col1_, col2_);ro.r(c)
    pval = float(rpy2.robjects.r('my.cor$p.value')[0])
    estimate = float(rpy2.robjects.r('my.cor$estimate')[0])
    return estimate, pval

def list_kw_dunn(names, lists, dep_var, dep_fac, outpath):
    dfs = []
    for item in range(len(names)):
        _ = pd.DataFrame()
        _[dep_var] = pd.Series(lists[item])
        _[dep_fac] = names[item]
        dfs.append(_)
    df1 = pd.concat(dfs).reset_index()
    del df1['index']

    #assert len(df1[df1['PEP'] < 0]) == 0

    data = KW_DUNN(df1, dep_var, dep_fac, outpath)


def mwu(dataframe, group1, group2):
    data = dataframe.copy()
    rdf =pandas2ri.py2ri(data)
    ro.globalenv['data'] = rdf
    c = 'print(data)';ro.r(c)
    c = "test = wilcox.test(data${} ~ data${})".format(group1,group2); ro.r(c)
    c = "print(test)";ro.r(c)
    c = "stat = as.list(test$statistic[1])[1]"; ro.r(c)
    #print
    c    = 'stat$W'
    stat = robjects.r(c)[0]

    c    ='test$p.value[1]'
    pval = robjects.r(c)[0]
    return stat, pval


def ANOVA_RM(dataframe, val, var, subject):   # pandas df, col to do anova on
    #c = "sink('{}')".format(outpath); ro.r(c)
    c = "library('FSA')"; ro.r(c)
    c = "library('nlme')"; ro.r(c)
    c = "require(multcomp)";ro.r(c)
    #var_dct = {}
    #count = 0
    #for v in list(set(dataframe[var])):
    #    var_dct[v] = str(count)
    #    count += 1
    #dataframe['tempvar'] = dataframe[var].apply(lambda x : var_dct[x])
    rdf =pandas2ri.py2ri(dataframe); ro.globalenv['data'] = rdf
    #c = "library(plyr)";ro.r(c)
    try:
        c = "your.bartlett = bartlett.test(data${}~data${})".format(val, var); ro.r(c)
        c = "print(your.bartlett)";ro.r(c)
        bartlett = robjects.r('your.bartlett$p.value')[0]
    except:
        bartlett = np.nan
    try:
        c = "your.fligner = fligner.test(data${}~data${})".format(val, var); ro.r(c)
        c = "print(your.fligner)";ro.r(c)
        fligner = robjects.r('your.fligner$p.value')[0]
    except:
        fligner = np.nan

    #print bartlett, fligner
    try:
        c = "your.lme =  lme({} ~ {}, data=data,random=~1|{}/{})".format(val, var, subject, var); ro.r(c)
        c = "your.anova = anova(your.lme)"; ro.r(c)
        c = "your.anova"
        an = pandas2ri.ri2py(ro.r[c])
        an = an.loc[var].to_dict()
        c=  "your.sum = summary(glht(your.lme, linfct=mcp({} = 'Tukey')), test = adjusted(type = 'bonferroni'))".format(var); ro.r(c)
        c =  "pvals = your.sum[10]$test$pvalues";ro.r(c)
        c = 'data.frame = data.frame(as.list(pvals))'; ro.r(c)
        c = 'data.frame'
        df = pandas2ri.ri2py(ro.r[c])
        df = df.reset_index()
        del df['index']
        df.loc[0, 'Tukey correction']          = 'Bonferroni'
        for key in an:
            df.loc[0, 'ANOVA '+key] = an[key]
    except:
        df = pd.DataFrame()

    return df

def KW_DUNN(dataframe, dep_var, dep_fac, outpath):   # pandas df, col to do anova on
    c = "sink('{}')".format(outpath); ro.r(c)
    c = "library('FSA')"; ro.r(c)
    rdf =pandas2ri.py2ri(dataframe); ro.globalenv['data'] = rdf
    c = "library(plyr)";ro.r(c)
    c = "your.kw = kruskal.test(data${} ~ data${})".format(dep_var, dep_fac); ro.r(c)
    kw_p = robjects.r('your.kw$p.value')[0]
    kw_statistic = robjects.r('your.kw$statistic')[0]
    kw_parameter = robjects.r('your.kw$parameter')[0]

    c = "your.dunn = dunnTest(data${}, data${}, kw=TRUE, two.sided = TRUE, method='bh')".format(dep_var,dep_fac);ro.r(c)
    c = "print(str(your.dunn))"; ro.r(c)
    dunn_correction = robjects.r('your.dunn$method')[0]
    c = "your.table = your.dunn$res"; ro.r(c)
    ro.r('print(your.table)')
    
    #df = pandas2ri.ri2py(ro.r[c])

    data = pd.DataFrame()
    data.loc[0, 'Dunn correction']       = dunn_correction
    data.loc[0, 'Kruskal-Wallis p-value'] = kw_p
    data.loc[0, 'Kruskal-Wallis statistic'] = kw_statistic
    data.loc[0, 'Kruskal-Wallis parameter'] = kw_parameter
    
    try:
        comp_dct = df.set_index('Comparison')['P.adj'].to_dict()
        for val in comp_dct:
            data.loc[0, val +' Dunn-test P.adj'] = comp_dct[val]
    except:
        pass
    
    return data

class up2ko:
    def __init__(self, up):
        self.up = up
        mg = mygene.MyGeneInfo()
        res = mg.query(self.up)
        self.ids = []
        self.names = []
        for i in res['hits']:
            symbol = i['entrezgene']
            ko = entrez2ko(symbol)
            self.ids.append(ko.ko)
            self.names.append(ko.name)
        self.ko = ';'.join(set(self.ids))
        self.name = ';'.join(set(self.names))

class string2ko: # uniprot id
    def __init__(self, string):
        c = 'library(KEGGREST)';ro.r(c)
        c = "res <- keggFind('ko', c('{}'))".format(string) ; ro.r(c)
        c = "print(res)" ; ro.r(c)
        c = 'names(res)'
        self.ko = ';'.join([i.split(':')[1] for i in robjects.r(c)])

        c = 'res'
        self.name= ';'.join(robjects.r(c))

class entrez2ko: # uniprot id
    def __init__(self, entrez):
        self.entrez = entrez
        c = 'library(KEGGREST)';ro.r(c)
        c = 'conv <- keggConv("genes", "ncbi-geneid:{}")'.format(self.entrez); ro.r(c) 
        c = 'ko <- keggGet(conv)'; ro.r(c)
        c = 'names(ko[[1]]$ORTHOLOGY)'
        self.ko = ';'.join(robjects.r(c))
        c = 'ko[[1]]$ORTHOLOGY'
        self.name= ';'.join(robjects.r(c))
        

    # keggConv("genes", "uniprot:Q05025"))
    # q <- keggGet("cjo:107318960")
    # q[[1]]$ORTHOLOGY   
    #c    = 'stat$W'
    #stat = robjects.r(c)[0]

def limma(df1, treated, control):   # pandas df, treated columns, control columns
    df = df1.copy()
    combined = treated + control 
    df[combined] =df[combined].astype(float)
    df  = df.reset_index()
    del df['index']
    names = {}
    ren = {}
    count = 1
    tr= []
    ct = []
    for i in df.columns:
        col = 'column{}'.format(count)
        names[col] = i
        ren[i] = col
        if i in treated:
            tr.append(col)
        elif i in control:
            ct.append(col)
        count += 1 
    df.rename(columns=ren, inplace = True)
    c="library(limma)"; ro.r(c)
    c="library(qvalue)"; ro.r(c)
    ctot = ct + tr
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c="str(data)"; #print ro.r(c)
    c= "data[ is.na(data) ] <- NA"; ro.r(c)
    c="tr <- c{}".format(tuple(tr)); ro.r(c); #print c
    c="ct <- c{}".format(tuple(ct)); ro.r(c); #print c
    c= "str(ct)"; #print ro.r(c)
    #c='''source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")'''
    #c='''source("source.functions.r.txt")'''
    ro.r(functions)
    c = 'data'
    design = []
    [design.append(2) for i in tr]
    [design.append(1) for i in ct]
    c="design <- model.matrix(~factor(c{}))".format(tuple(design))
    #print c
    ro.r(c)
    #print ro.r('str(data)')
    c='''colnames(design) <- c("Intercept", "Diff")'''
    ro.r(c)
    c='res.eb <- eb.fit(data[, c(tr,ct)], design)'
    ro.r(c)
    c="res.eb"
    #print ro."r(c)
    dfebi = pandas2ri.ri2py(ro.r[c])
    assert 'index' not in dfebi.columns
    dfebi.index = dfebi.reset_index()['index'].apply(int)
    df = df.join(dfebi)
    df = df.rename(columns=names)
    df = df.sort('p.mod',ascending = True)
    return df

def trt(df1,treated, control):   # pandas df, treated columns, control columns
    df = df1.copy()
    combined = treated + control
    df[combined] =df[combined].astype(float)
    df  = df.reset_index()
    del df['index']
    names = {}
    ren = {}
    count = 1
    tr= []
    ct = []
    for i in df.columns:
        col = 'column{}'.format(count)
        names[col] = i
        ren[i] = col
        if i in treated:
            tr.append(col)
        elif i in control:
            ct.append(col)
        count += 1
    df.rename(columns=ren, inplace = True)
    c="library(limma)"; ro.r(c)
    c="library(qvalue)"; ro.r(c)
    ctot = ct + tr
    rdf =pandas2ri.py2ri(df)
    ro.globalenv['data'] = rdf
    c="str(data)"; #print ro.r(c)
    c= "data[ is.na(data) ] <- NA"; ro.r(c)
    c="tr <- c{}".format(tuple(tr)); ro.r(c); #print c
    c="ct <- c{}".format(tuple(ct)); ro.r(c); #print c
    c= "str(ct)"; #print ro.r(c)
    #c='''source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")'''
    c='''source("source.functions.r.txt")'''
    ro.r(c)
    c = 'data'
    design = []
    [design.append(2) for i in tr]
    [design.append(1) for i in ct]
    c="design <- model.matrix(~factor(c{}))".format(tuple(design))
    #print c
    ro.r(c)
    #print ro.r('str(data)')
    c='''colnames(design) <- c("Intercept", "Diff")'''
    ro.r(c)
    c='res.eb <- trt.fit(data[, c(tr,ct)], design)'
    ro.r(c)
    c="res.eb"
    #print ro."r(c)
    dfebi = pandas2ri.ri2py(ro.r[c])
    assert 'index' not in dfebi.columns
    dfebi.index = dfebi.reset_index()['index'].apply(int)
    df = df.join(dfebi)
    df = df.rename(columns=names)
    df = df[df['adj.P.Val'].notnull()]
    df = df.sort('adj.P.Val',ascending = True)
    return df
