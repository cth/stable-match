# Experiments with GSA matching 

source("match.r")
library(VennDiagram)
library(grid)

sample.case.control.pop <- function(nfeatures,ncases,ncontrols) {
	sample.pop <- data.frame(id=1:(ncases+ncontrols))

	for(i in 1:nfeatures)
		sample.pop[[i]] = rnorm(ncases + ncontrols)

	list(
		 cases = sample.pop[1:ncases,],
		 controls = sample.pop[seq(ncases+1,ncontrols+ncases),],
		 match.features = names(sample.pop)[2:nfeatures]
	)
}
measurements <- c("longest.distance", "mean.distance","median.distance", "sd.distance")


run.increasing.cmp <- function(num.experiments=10,num.features=2,control.population=50,num.cases.seq=seq(5,50,5)) { 
	runs <- list() 

	for(i in 1:num.experiments) {
		results <- list()
		results$gsa <- list()
		results$opt <- list()
		results$nn <- list()
		results$measure.points <- c()

		for(ncases in num.cases.seq) {
			results$measure.points <- c(results$measure.points,ncases)
			pop <- sample.case.control.pop(num.features,ncases,control.population) 
			for(algo in c("gsa", "opt","nn")) {
				match <- match.cc(pop$cases,pop$controls,pop$match.features,id.column=1,distance.measure="mahalanobis",match.method=algo)

				for(m in measurements) 
					results[[algo]][[m]] <- c( results[[algo]][[m]], match$test.metrics[[m]])
			}
		}
		runs[[length(runs)+1]] <- results
	}

	plot.measurement <- function(name, gsa, lp, nn, measure.points,ylab="") {
		pdf(paste0("plots/", gsub("\\.", "_", name),".pdf"))
		plot(measure.points, gsa, ylim=range(gsa,lp,nn),col="blue",type="b",main=gsub("\\.", " ", name), xlab="#cases",ylab=ylab) 
		lines(measure.points,lp,col="red",t="b")
		lines(measure.points,nn,col="green",t="b")
		legend("topleft", c("gsa", "opt","nn"), col=c("blue", "red","green"), pch=21)
		dev.off()
	}

	matplot.measurement <- function(name, gsa, lp,nn, measure.points,ylab="") {
		pdf(paste0("plots", "mat.", name,".pdf"))
		print(dim(rbind(gsa,lp)))
		color=c(rep("blue",nrow(gsa)),rep("red", nrow(lp)), rep("green", nrow(nn)))
		matplot(measure.points, t(rbind(gsa,lp,nn)), ylim=range(gsa,lp,nn),col=color,type="b",main=name, xlab="#cases",ylab=ylab) 
		legend("topleft", c("gsa", "opt", "nn"), col=c("blue", "red","green"), pch=21)
		dev.off()
	}


	for(m in measurements) {
		gsa <- matrix(unlist(lapply(runs, function(run) { run[["gsa"]][[m]] })),length(runs),byrow=T)
		lp <- matrix(unlist(lapply(runs, function(run) { run[["opt"]][[m]] })),length(runs),byrow=T)
		nn <- matrix(unlist(lapply(runs, function(run) { run[["nn"]][[m]] })),length(runs),byrow=T)
		gsa.means <- apply(gsa,2,mean)
		lp.means <- apply(lp,2,mean)
		nn.means <- apply(nn,2,mean)
		plot.measurement(m, gsa.means, lp.means, nn.means, results$measure.points)
#		matplot.measurement(m, gsa, lp, nn, results$measure.points)
	}
}

run.cmp.one <- function(num.experiments=100,num.cases=10,control.population=30,comparison.metric="mean") {
	runs <- list() 

	venn.area <- list()
	distances <- list()


	for(algo in c("gsa", "opt", "nn")) {
		distances[[algo]] <- c()	
	}

	for(m in measurements) {
		venn.area[[m]] <- list()
		venn.area[[m]]$area1 <- 0
		venn.area[[m]]$area2 <- 0
		venn.area[[m]]$area3 <- 0
		venn.area[[m]]$n12 <- 0
		venn.area[[m]]$n23 <- 0
		venn.area[[m]]$n13 <- 0
		venn.area[[m]]$n123 <- 0
	}

	for(i in 1:num.experiments) {
		results <- list()
		results$gsa <- list()
		results$opt <- list()
		results$measure.points <- c()

		pop <- sample.case.control.pop(2,num.cases,control.population) 
		for(algo in c("gsa", "opt", "nn")) {
			results[[algo]] <- list()
			match <- match.cc(pop$cases,pop$controls,pop$match.features,id.column=1,distance.measure="mahalanobis",match.method=algo)
			distances[[algo]] <- c(distances[[algo]],match$distances$control1)
			for(m in measurements) 
				results[[algo]][[m]] <- match$test.metrics[[m]]
		}

		for(m in measurements) {
			min <- 1000
			best <- c()
			for (method in c("gsa", "opt", "nn")) {
				if ( results[[method]][[m]] < min ) {
					min <- results[[method]][[m]] 
					best <- method
				} else if( results[[method]][[m]] == min) {
					best <- c(best,method)
				}
			}

			venn.area[[m]][["area1"]] <- venn.area[[m]][["area1"]] + sum("gsa" %in% best)
			venn.area[[m]][["area2"]] <- venn.area[[m]][["area2"]] + sum("opt" %in% best)
			venn.area[[m]][["area3"]] <- venn.area[[m]][["area3"]] + sum("nn" %in% best)
			venn.area[[m]][["n12"]] <- venn.area[[m]][["n12"]] + sum("gsa" %in% best & "opt" %in% best)
			venn.area[[m]][["n23"]] <- venn.area[[m]][["n23"]] + sum("opt" %in% best & "nn" %in% best)
			venn.area[[m]][["n13"]] <- venn.area[[m]][["n13"]] + sum("gsa" %in% best & "nn" %in% best)
			venn.area[[m]][["n123"]] <- venn.area[[m]][["n123"]] + sum("gsa" %in% best & "opt" %in% best & "nn" %in% best)
		} 
	}

	pdf(paste0("plots/","venn_",gsub("\\.","_",m),".pdf")) 
	venn.plot <- draw.triple.venn(
		area1 = venn.area[[m]][["area1"]],
		area2 = venn.area[[m]][["area2"]],
		area3 = venn.area[[m]][["area3"]],
		n12 = venn.area[[m]][["n12"]],
		n23 = venn.area[[m]][["n23"]],
		n13 = venn.area[[m]][["n13"]],
		n123 = venn.area[[m]][["n123"]],
		category = c("gsa", "opt", "nn"),
		cex=2,
		#fill = c("blue", "red", "green"),
		#cat.col = c("blue", "red", "green"),
		scaled = T)
	grid.draw(venn.plot)
	dev.off()
	grid.newpage()

	print("------------------------------------------------------")
	print(distances)
	pdf("plots/qqplot_gsa_opt.pdf")
	plot(sort(distances[["gsa"]]),sort(distances[["opt"]]))
	abline(0,1,col="red")
	dev.off()

	pdf("plots/qqplot_gsa_nn.pdf")
	plot(sort(distances[["gsa"]]),sort(distances[["nn"]]))
	abline(0,1,col="red")
	dev.off()
	
	pdf("plots/qqplot_opt_nn.pdf")
	plot(sort(distances[["opt"]]),sort(distances[["nn"]]))
	abline(0,1,col="red")
	dev.off()
}

run.cmp.one(num.experiments=100)
#run.increasing.cmp()
