# Experiments with GSA matching 

source("match.r")

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

runs <- list() 
num.experiments <- 5 

for(i in 1:num.experiments) {
	results <- list()
	results$gsa <- list()
	results$lpsolve <- list()
	results$measure.points <- c()

	for(ncases in seq(5,50,5)) {
		results$measure.points <- c(results$measure.points,ncases)
		for(algo in c("gsa", "lpsolve")) {
			pop <- sample.case.control.pop(2,ncases,50) 

			match <- match.cc(pop$cases,pop$controls,pop$match.features,id.column=1,distance.measure="mahalanobis",match.algorithm=algo)

			for(m in measurements) 
				results[[algo]][[m]] <- c( results[[algo]][[m]], match$test.metrics[[m]])
		}
	}
	runs[[length(runs)+1]] <- results
}

plot.measurement <- function(name, gsa, lp, measure.points,ylab="") {
	pdf(paste0(name,".pdf"))
	plot(measure.points, gsa, ylim=range(gsa,lp),col="blue",type="b",main=name, xlab="#cases",ylab=ylab) 
	lines(measure.points,lp,col="red",t="b")
	legend("topleft", c("gsa", "lp"), col=c("blue", "red"), pch=21)
	dev.off()
}

matplot.measurement <- function(name, gsa, lp, measure.points,ylab="") {
	pdf(paste0("mat.", name,".pdf"))
	print(dim(rbind(gsa,lp)))
	color=c(rep("blue",nrow(gsa)),rep("red", nrow(lp)))
	matplot(measure.points, t(rbind(gsa,lp)), ylim=range(gsa,lp),col=color,type="b",main=name, xlab="#cases",ylab=ylab) 
	legend("topleft", c("gsa", "lp"), col=c("blue", "red"), pch=21)
	dev.off()
}


#print(lapply(runs, function(run) { run[["gsa"]][["longest.distance"]] }))
#tmp <- matrix(unlist(lapply(runs, function(run) { run[["gsa"]][["longest.distance"]] })),length(runs),byrow=T)
#print("----------------------")
#print(tmp)
#means <- apply(tmp,2,mean)


for(m in measurements) {
	gsa <- matrix(unlist(lapply(runs, function(run) { run[["gsa"]][["longest.distance"]] })),length(runs),byrow=T)
	lp <- matrix(unlist(lapply(runs, function(run) { run[["lpsolve"]][["longest.distance"]] })),length(runs),byrow=T)
	gsa.means <- apply(gsa,2,mean)
	lp.means <- apply(lp,2,mean)
	plot.measurement(m, gsa.means, lp.means,results$measure.points)
	matplot.measurement(m, gsa, lp,results$measure.points)
}
