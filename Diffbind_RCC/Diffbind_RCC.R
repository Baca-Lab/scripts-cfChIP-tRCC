library(DiffBind) #v3.10.1


#generate dba objects
DiffRCC<- dba(sampleSheet = "./manifestDiffbindRCC_example.csv")

DiffRCCcounts<- dba.count(DiffRCC)
DiffRCCAnalyse<-dba.analyze(DiffRCCcounts)

#generate report
DiffRCC_report <- dba.report(DiffRCCAnalyse)

#generate PCA and Heatmap
dba.plotPCA(DiffRCCcounts,attributes=DBA_CONDITION,
            label=DBA_ID,
            labelSize=0.5, 
            vColors=c("#B09C8599","#DC000099"),
            dotSize=1)

dba.plotHeatmap(DiffRCCcounts,
                DBA_TISSUE,
                label=DBA_CONDITION,
                colSideCols=c("#B09C8599","#DC000099"))






