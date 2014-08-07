decompose.formula <- function(formula,bindings=list(),lhs=T) {
	if (class(formula)=="formula" | class(formula)=="call") {
		lf <- as.list(formula)
		if (length(lf)==3)
			c(bindings,decompose.formula(lf[[2]],list(),lhs), decompose.formula(lf[[3]],list(),ifelse(lf[[1]]=="|",lhs,F)))
		else
			c(bindings,decompose.formula(lf[[2]],list(),lhs))
	} else if (class(formula)=="name") {
			bindings[[as.character(formula)]]=ifelse(lhs,"mvar","stvar")
			bindings
	} else {
		NULL
	}
}


b <- decompose.formula(~ bmi + age | factor) 

