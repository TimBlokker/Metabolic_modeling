library(lattice)
data<-read.csv('plotfile.csv')
attach(data)
fluxes <- data.frame(metabolite, flux_mmcb, flux_mgam, flux_rcm)
lbs <- data.frame(metabolite, lower_bound_mmcb, lower_bound_mgam, lower_bound_rcm)
nonzero.fluxes <- fluxes[as.logical(rowSums(fluxes[,2:4] != 0)), ]
nonzero.lbs <- lbs[as.logical(rowSums(lbs[,2:4] != 0)), ]

full.fluxes <- nonzero.fluxes[nonzero.fluxes$flux_mmcb != 0 & nonzero.fluxes$flux_mgam != 0 & nonzero.fluxes$flux_rcm != 0,]

xyplot(factor(nonzero.fluxes$metabolite)~nonzero.fluxes$flux_mmcb+nonzero.fluxes$flux_mgam+nonzero.fluxes$flux_rcm, 
        main="Flux/Metabolite", 
        ylab="Flux", 
        xlab="Metabolite",
        col=c("red","green","blue"),
        xlim=c(-1000,1000),
        scales=list(y=list(cex=0.4)))

xyplot(factor(nonzero.lbs$metabolite)~nonzero.lbs$lower_bound_mmcbb+nonzero.lbs$lower_bound_mgam+nonzero.lbs$lower_bound_rcm, 
       main="Flux/Metabolite", 
       ylab="Flux", 
       xlab="Metabolite",
       col=c("red","green","blue"),
       xlim=c(-1000,1000),
       scales=list(y=list(cex=0.4, tck=2)))

xyplot(factor(full.fluxes$metabolite)~full.fluxes$flux_mmcb+full.fluxes$flux_mgam+full.fluxes$flux_rcm, 
        main="Flux per Metabolite", 
        ylab="Metabolite", 
        xlab="Flux",
        col=c("red","green","blue"),
        pch=20,
        scales=list(y=list(cex=0.4), x=list(limits=c(-1000,1000))),
        auto.key = list(points=F, col=c("red","green","blue"), 
                                text = c("mMCB", "mGAM", "RCM"), title = "Medium", 
                                x=0.2, y=0.9))

xyplot(factor(nonzero.fluxes$metabolite)~nonzero.fluxes$flux_rcm+nonzero.fluxes$flux_mgam+nonzero.fluxes$flux_mmcb, 
                main="Flux per Metabolite", 
                ylab="Metabolite", 
                xlab="Flux",
                col=c("red","green","blue"),
                pch=c(15,16,17),
                xlim=c(-1000,1000),
                scales=list(y=list(cex=0.4)),
                auto.key = list(points=F, col=c("red","green","blue"), 
                                                        text = c("RCM", "mGAM", "mMCB"), title = "Medium",
                                                        x=0.2, y=0.8))

xyplot(factor(nonzero.lbs$metabolite)~nonzero.lbs$lower_bound_rcm+nonzero.lbs$lower_bound_mgam+nonzero.lbs$lower_bound_mmcb, 
       main="Lower Bound per Metabolite", 
       ylab="Metabolite", 
       xlab="Lower Bound",
       col=c("red","green","blue"),
       pch=c(15,16,17),
       xlim=c(-1000,0),
       scales=list(y=list(cex=0.4)),
       auto.key = list(points=F, col=c("red","green","blue"), 
                       text = c("RCM", "mGAM", "mMCB"), title = "Medium",
                       x=0.1, y=0.9))
