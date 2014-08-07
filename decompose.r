# radical idea:
# xpath like expressions for R (Rpath)

decompose.formula <- function(formula,bindings=list(),lhs=T) {
	print(paste("decompose.formula(", as.character(formula), ",", class(formula),")"))
	if (class(formula)=="formula" | class(formula)=="call") {
		#print(">>")
		lf <- as.list(formula)
		#print(lf)
		# might be able to simplify:
		if (lf[[1]] == "|") {
			if (length(lf)==3)
				c(bindings,decompose.formula(lf[[2]],list(),lhs), decompose.formula(lf[[3]],list(),F))
			else
				c(bindings,decompose.formula(lf[[2]],list(),lhs))
		} else {
			if (length(lf)==3)
				c(bindings,decompose.formula(lf[[2]],list(),lhs), decompose.formula(lf[[3]],list(),lhs))
			else
				c(bindings,decompose.formula(lf[[2]],list(),lhs))
		}
	} else if (class(formula)=="name") {
			print(as.character(formula))
			bindings[[as.character(formula)]]=ifelse(lhs,"mvar","stvar")
			print(bindings)
	} else {
		NULL
	}
}


b <- decompose.formula(~ bmi + age | factor) 

