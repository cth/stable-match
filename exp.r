# Experiments with GSA matching 

source("match.r")

sample.case.control.pop <- function(nfeatures,ncases,ncontrols) {
	sample.pop <- data.frame(id=1:(ncases+ncontrols))

	for(i in 1:nfeatures)
		sample.pop[[i]] = rnorm(ncases + ncontrols) 
#		sample.pop[[i]] = rnorm(ncases + ncontrols, 
#								mean=sample(0:nfeatures,1),
#								sd=sample(1:nfeatures,1))

	list(
		 cases = sample.pop[1:ncases,],
		 controls = sample.pop[seq(ncases+1,ncontrols+ncases),],
		 match.features = names(sample.pop)[2:nfeatures]
	)
}

times.gsa <- c()
times.lp <- c()
longest.gsa <- c()
longest.lp <- c()
mean.gsa <- c()
mean.lp <- c()
median.gsa <- c()
median.lp <- c()
measure.points <- c()

results <- list()
results$gsa <- list()
results$lpsolve <- list()

for(ncases in seq(2,20,2)) {
	measure.points <- c(measure.points,ncases)
	for(algo in c("gsa", "lpsolve")) {
		pop <- sample.case.control.pop(3,ncases,20) 

		match <- match.cc(pop$cases,pop$controls,pop$match.features,id.column=1,distance.measure="mahalanobis",match.algorithm=algo)

		
		if (is.null(results[[algo]]$times)
			results[[algo]]$times <- c()
		else
			results[[algo]]$times <- c( results[[algo]]$times, match$test.metrics$elapsed.time[3]))

		if (is.null(results[[algo]]$longest)
			results[[algo]]$longest <- c()
		else
			results[[algo]]$longest <- c( results[[algo]]$longest, match$test.metrics$longest.distance))

		if (is.null(results[[algo]]$mean)
			results[[algo]]$mean <- c()
		else
			results[[algo]]$mean <- c( results[[algo]]$median, match$test.metrics$mean.distance))

		if (is.null(results[[algo]]$median)
			results[[algo]]$median <- c()
		else
			results[[algo]]$median <- c( results[[algo]]$median, match$test.metrics$median.distance))


			
		results[[algo]]$times <- 
		times.gsa = c(times.gsa, match.gsa$test.metrics$elapsed.time[3])

		match.gsa <- match.cc(pop$cases,pop$controls,pop$match.features,id.column=1,distance.measure="mahalanobis",match.algorithm="gsa")
		times.gsa = c(times.gsa, match.gsa$test.metrics$elapsed.time[3])
		longest.gsa = c(longest.gsa, match.gsa$test.metrics$max.distance)
		mean.gsa = c(mean.gsa, match.gsa$test.metrics$mean.distance)
		median.gsa = c(median.gsa, match.gsa$test.metrics$median.distance) 

		match.lp <- match.cc(pop$cases,pop$controls,pop$match.features,id.column=1,distance.measure="mahalanobis",match.algorithm="lpsolve")
		times.lp = c(times.lp, match.lp$test.metrics$elapsed.time[3])
		longest.lp = c(longest.lp, match.lp$test.metrics$max.distance)
		mean.lp = c(mean.lp, match.lp$test.metrics$mean.distance)
		median.lp = c(median.lp, match.lp$test.metrics$median.distance) 
	}
}

plot.measurement <- function(name, gsa, lp, y.points) {
	print(paste("plot", name))

	min.x <- min(c(gsa,lp))
	max.x <- max(c(gsa,lp))
	print(min.x)
	print(max.x)
	pdf(paste0(name,".pdf"))
	plot(gsa, y.points, ylim=c(min.x,max.x),col="blue",type="b") 
	lines(lp,y.points,col="red")
	legend("topleft", c("gsa", "lp"), col=c("blue", "red"), pch=21)
	dev.off()
}

print(measure.points)

plot.measurement("times", times.gsa,times.lp,measure.points)
plot.measurement("longest", times.gsa,times.lp,measure.points)
#plot("longest", longest.gsa,times.lp)
#plot("mean", mean.gsa,times.lp)
#plot("median", median.gsa,times.lp)

#pdf("times.pdf")
#plot(times.gsa,col="blue", type="l")
#lines(times.lp,col="red")
#dev.off()
#
#pdf("longest.pdf")
#plot(longest.gsa,col="blue", type="l")
#lines(longest.lp,col="red")
#dev.off()
#
#pdf("mean.pdf")
#plot(mean.gsa,col="blue", type="l")
#lines(mean.lp,col="red")
#dev.off()
#
#pdf("median.pdf")
#plot(median.gsa,col="blue", type="l")
#lines(median.lp,col="red")
#dev.off()
