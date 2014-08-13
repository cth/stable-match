# Experiments with GSA matching 

source("match.r")

sample.case.control.pop <- function(nfeatures,ncases,ncontrols) {
	sample.pop <- data.frame(id=1:(ncases+ncontrols))

	for(i in 1:nfeatures)
		sample.pop[[i]] = sin(exp(rnorm(ncases + ncontrols)))

	list(
		 cases = sample.pop[1:ncases,],
		 controls = sample.pop[seq(ncases+1,ncontrols+ncases),],
		 match.features = names(sample.pop)[2:nfeatures]
	)
}
measurements <- c("longest.distance", "mean.distance","median.distance", "sd.distance")


run.increasing.cmp <- function(num.experiments=5,num.features=2,control.population=50,num.cases.seq=seq(5,50,5)) { 
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
		pdf(paste0("plots/", name,".pdf"))
		plot(measure.points, gsa, ylim=range(gsa,lp,nn),col="blue",type="b",main=name, xlab="#cases",ylab=ylab) 
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
		matplot.measurement(m, gsa, lp, nn, results$measure.points)
	}
}

run.cmp.one <- function(num.experiments=20,num.cases=20,control.population=60,comparison.metric="mean") {
	runs <- list() 

	wins <- list()
	for(m in measurements) {
		wins[[m]] <- 0
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
			for(m in measurements) 
				results[[algo]][[m]] <- match$test.metrics[[m]]
		}

		for(m in measurements) {
			print("----------------------------")
			print(results[["gsa"]])
			if (results[["gsa"]][[m]] < results[["opt"]][[m]]) {
				print("gsa wins")
				wins[[m]] <- wins[[m]] + 1
			} else {
				print("lpsolve wins")
			}
		}
	} 

	print("##################### Score: ############################")
	print(wins)
}

#run.cmp.one()


run.increasing.cmp()
