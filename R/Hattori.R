###############################################################################
#
# Probabilistic Representation Theory for syllogistic reasoning
#
#   Version 4.0
#
#   (c)2016 Masasi Hattori
#
###############################################################################
#
# Usage:
#
#  PRT(data, n=7)       Probabilistic Representation Theory
#  PHM(data)            Probability Heuristic Model
#  pMM(data)            Parameterized model based on Mental Medel Theory
#
# Arguments:
#
#  data   Dataframe for experimental data
#  n      The number of sample in Sample Mental Model
#
# Details:
#
#  The format of "data" is expected to be 64 x 5 as follows:
#
#               A           I           E           O          N
#  AA1 0.878048780 0.046747967 0.008130081 0.000000000 0.06707317
#  AA2 0.563265306 0.061224490 0.006122449 0.002040816 0.36734694
#  ...
#
# Examples:
#
#  > data1 <- read.csv("C:/R-files/Hypothesis/HattoriData/CO99.csv", header=T, row.names=1)
#
#  > prt <- PRT(data)
#  > model <- prt$m
#  > x  <- prt$p[1]; c  <- prt$p[2]
#  > ua <- prt$p[3]; ui <- prt$p[4]; ue <- prt$p[5]; uo <- prt$p[6]
#
#  > rmsd <- RMSD.syl(data, model)
#  > corr <- CORR.syl(data, model)
#  > mean(rmsd)
#  > fisher_inv(mean(fisher(corr)))
#
#  > phm <- PHM(data)
#  > model <- phm$m
#  > pa <- phm$p[1]; pi <- phm$p[2]; pe <- phm$p[3]; po <- phm$p[4]
#  > pent <- phm$p[5]; perr <- phm$p[6]
#
#  > pmm <- pMM(data)
#  > model <- pmm$m
#  > p1m1 <- par[1]; p2m1 <- par[2]; p2m2 <- par[3]
#  > p3m1 <- par[4]; p3m2 <- par[5]; p3m3 <- par[6]
#
#------------------------------------------------------------------------------
# Usage:
#
#  PRT.AMFO(data, n=7)
#  PRT.MFIE(data, n=7)
#
# Arguments:
#
#  data   Dataframe for experimental data
#  n      The number of sample in Sample Mental Model
#
# Details:
#
#  The format of "data" for PRT.AMFO is expected to be 64 x 5 as follows:
#
#         A    M    F    O    N
#  AA1 0.85 0.00 0.00 0.10 0.05
#  AA2 0.60 0.05 0.10 0.05 0.20
#  ...
#  AM1 0.05 0.85 0.00 0.10 0.00
#  ...
#
#  The format of "data" for PRT.MFIE is expected to be 64 x 5 as follows:
#
#         M    F    I    E    N
#  MM1 0.45 0.10 0.20 0.15 0.20
#  MM2 0.50 0.15 0.05 0.05 0.25
#  ...
###############################################################################

MD <- 
  c("AA","AI","IA","AE","EA","AO","OA",
    "II","IE","EI","IO","OI",
    "EE","EO","OE",
    "OO")

MD.AMFO <- 
  c("AA","AM","MA","AF","FA","AO","OA",
    "MM","MF","FM","MO","OM",
    "FF","FO","OF",
    "OO")

MD.MFIE <- 
  c("MM","MF","FM","MI","IM","ME","EM",
    "FF","FI","IF","FE","EF",
    "II","IE","EI",
    "EE")

###############################################################################

#
# To return a predicted probability distribution of responses
#
Syl <- function(mod, fig, x, c, n, ua, ui, ue, uo)
{
	#################################
	# STEP 1 & 2 OF THE MODEL

	AMFIEO <- Syl.simple(mod, fig, x, c, n)
	AIEO <- AMFIEO[-c(2,3)]

	#################################
	# STEP 3 OF THE MODEL

	# Probabilities to pass a single n-th test
	# PrX1 = 1st test, PrX2 = 2nd test
	# i.e., 0 <= (PrX1, PrX2, PrX3, PrX4) <= 1
	# The min-heuristic:
	# The least informative premise defines the order of the test
	if (mod == "AA") {
		# min = A
		PrX1 <- AIEO[1]
		PrX2 <- AIEO[2]
		PrX3 <- AIEO[3]
		PrX4 <- AIEO[4]
	} else if (mod == "AI" || mod == "IA" || mod == "II") {
		# min = I
		PrX1 <- AIEO[2]
		PrX2 <- AIEO[3]
		PrX3 <- AIEO[4]
		PrX4 <- AIEO[1]
	} else if (mod == "AE" || mod == "EA" ||
		 mod == "IE" || mod == "EI" || mod == "EE") {
		# min = E
		PrX1 <- AIEO[3]
		PrX2 <- AIEO[4]
		PrX3 <- AIEO[1]
		PrX4 <- AIEO[2]
	} else if (mod == "AO" || mod == "OA" ||
		 mod == "IO" || mod == "OI" ||
		 mod == "EO" || mod == "OE" || mod == "OO") {
		# min = O
		PrX1 <- AIEO[4]
		PrX2 <- AIEO[1]
		PrX3 <- AIEO[2]
		PrX4 <- AIEO[3]
	}

	# Probabilities that the conclusion is validated by the n-th test
	# PX1 = 1st test, PX2 = 2nd test, ...
	# i.e., PX1 + PX2 + PX3 + PX4 + residual = 1
	PX1 <- PrX1					# pass
	PX2 <- (1 - PrX1)*PrX2				# fail -> pass
	PX3 <- (1 - PrX1)*(1 - PrX2)*PrX3		# fail*2 -> pass
	PX4 <- (1 - PrX1)*(1 - PrX2)*(1 - PrX3)*PrX4	# fail*3 -> pass

	if (mod == "AA") {
		# min = A
		# A <- PX1; I <- PX2; E <- PX3; O <- PX4
		AIEO <- c(PX1, PX2, PX3, PX4)
	} else if (mod == "AI" || mod == "IA" || mod == "II") {
		# min = I
		# I <- PX1; E <- PX2; O <- PX3; A <- PX4
		AIEO <- c(PX4, PX1, PX2, PX3)
	} else if (mod == "AE" || mod == "EA" ||
		 mod == "IE" || mod == "EI" || mod == "EE") {
		# min = E
		# E <- PX1; O <- PX2; A <- PX3; I <- PX4
		AIEO <- c(PX3, PX4, PX1, PX2)
	} else if (mod == "AO" || mod == "OA" ||
		 mod == "IO" || mod == "OI" ||
		 mod == "EO" || mod == "OE" || mod == "OO") {
		# min = O
		# O <- PX1; A <- PX2; I <- PX3; E <- PX4
		AIEO <- c(PX2, PX3, PX4, PX1)
	}

	# The max-heuristic:
	# The most informative premise determines the confidence in
	# the conclusion
	if (mod == "AA" || 
	    mod == "AI" || mod == "IA" ||
	    mod == "AE" || mod == "EA" ||
	    mod == "AO" || mod == "OA") {
		# max = A
		AIEO <- AIEO*ua
	} else if (mod == "II" ||
		 mod == "IE" || mod == "EI" ||
		 mod == "IO" || mod == "OI") {
		# max = I
		AIEO <- AIEO*ui
	} else if (mod == "EE" ||
		 mod == "EO" || mod == "OE") {
		# max = E
		AIEO <- AIEO*ue
	} else if (mod == "OO") {
		# max = O
		AIEO <- AIEO*uo
	}

	return(c(AIEO, 1 - sum(AIEO)))
}


Syl.AMFO <- function(mod, fig, x, c, n, ua, um, uf, uo, m, f)
{
	#################################
	# STEP 1 & 2 OF THE MODEL

	AMFIEO <- Syl.simple(mod, fig, x, c, n, m, f)
	AMFO <- AMFIEO[-c(4,5)]

	#################################
	# STEP 3 OF THE MODEL

	# Probabilities to pass a single n-th test
	# PrX1 = 1st test, PrX2 = 2nd test
	# i.e., 0 <= (PrX1, PrX2, PrX3, PrX4, PrX5, PrX6) <= 1
	# The min-heuristic:
	# The least informative premise defines the order of the test

	if (mod == "AA") {
		# min = A
		PrX1 <- AMFO[1]
		PrX2 <- AMFO[2]
		PrX3 <- AMFO[3]
		PrX4 <- AMFO[4]
	} else if (mod == "AM" || mod == "MA" || mod == "MM") {
		# min = M
		PrX1 <- AMFO[2]
		PrX2 <- AMFO[3]
		PrX3 <- AMFO[4]
		PrX4 <- AMFO[1]
	} else if (mod == "AF" || mod == "FA" ||
		   mod == "MF" || mod == "FM" || mod == "FF") {
		# min = F
		PrX1 <- AMFO[3]
		PrX2 <- AMFO[4]
		PrX3 <- AMFO[1]
		PrX4 <- AMFO[2]
	} else if (mod == "AO" || mod == "OA" || mod == "MO" || mod == "OM" || 
		   mod == "FO" || mod == "OF" || mod == "OO") {
		# min = O
		PrX1 <- AMFO[4]
		PrX2 <- AMFO[1]
		PrX3 <- AMFO[2]
		PrX4 <- AMFO[3]
	}

	# Probabilities that the conclusion is validated by the n-th test
	# PX1 = 1st test, PX2 = 2nd test, ...
	# i.e., PX1 + PX2 + PX3 + PX4 + residual = 1
	PX1 <- PrX1					# pass
	PX2 <- (1 - PrX1)*PrX2				# fail -> pass
	PX3 <- (1 - PrX1)*(1 - PrX2)*PrX3		# fail*2 -> pass
	PX4 <- (1 - PrX1)*(1 - PrX2)*(1 - PrX3)*PrX4	# fail*3 -> pass

	if (mod == "AA") {
		# min = A
		# A <- PX1; M <- PX2; F <- PX3; O <- PX4
		AMFO <- c(PX1, PX2, PX3, PX4)
	} else if (mod == "AM" || mod == "MA" || mod == "MM") {
		# min = M
		# M <- PX1; F <- PX2; O <- PX3; A <- PX4
		AMFO <- c(PX4, PX1, PX2, PX3)
	} else if (mod == "AF" || mod == "FA" ||
		   mod == "MF" || mod == "FM" || mod == "FF") {
		# min = F
		# F <- PX1; O <- PX2; A <- PX3; M <- PX4
		AMFO <- c(PX3, PX4, PX1, PX2)
	} else if (mod == "AO" || mod == "OA" ||
		 mod == "MO" || mod == "OM" ||
		 mod == "FO" || mod == "OF" || mod == "OO") {
		# min = O
		# O <- PX1; A <- PX2; M <- PX3; F <- PX4
		AMFO <- c(PX2, PX3, PX4, PX1)
	}

	# The max-heuristic:
	# the most informative premise determines the confidence in
	# the conclusion
	if (mod == "AA" || 
	    mod == "AM" || mod == "MA" ||
	    mod == "AF" || mod == "FA" ||
	    mod == "AO" || mod == "OA") {
		# max = A
		AMFO <- AMFO*ua
	} else if (mod == "MM" ||
		 mod == "MF" || mod == "FM" ||
		 mod == "MO" || mod == "OM") {
		# max = M
		AMFO <- AMFO*um
	} else if (mod == "FF" ||
		 mod == "FO" || mod == "OF") {
		# max = F
		AMFO <- AMFO*uf
	} else if (mod == "OO") {
		# max = O
		AMFO <- AMFO*uo
	}

	return(c(AMFO, 1 - sum(AMFO)))
}


Syl.MFIE <- function(mod, fig, x, c, n, um, uf, ui, ue, m, f)
{
	#################################
	# STEP 1 & 2 OF THE MODEL

	AMFIEO <- Syl.simple(mod, fig, x, c, n, m, f)
	MFIE <- AMFIEO[-c(1,6)]

	#################################
	# STEP 3 OF THE MODEL

	# Probabilities to pass a single n-th test
	# PrX1 = 1st test, PrX2 = 2nd test
	# i.e., 0 <= (PrX1, PrX2, PrX3, PrX4, PrX5, PrX6) <= 1
	# The min-heuristic:
	# The least informative premise defines the order of the test

	if (mod == "MM") {
		# min = M
		PrX1 <- MFIE[1]
		PrX2 <- MFIE[2]
		PrX3 <- MFIE[3]
		PrX4 <- MFIE[4]
	} else if (mod == "MF" || mod == "FM" || mod == "FF") {
		# min = F
		PrX1 <- MFIE[2]
		PrX2 <- MFIE[3]
		PrX3 <- MFIE[4]
		PrX4 <- MFIE[1]
	} else if (mod == "MI" || mod == "IM" ||
		   mod == "FI" || mod == "IF" || mod == "II") {
		# min = I
		PrX1 <- MFIE[3]
		PrX2 <- MFIE[4]
		PrX3 <- MFIE[1]
		PrX4 <- MFIE[2]
	} else if (mod == "ME" || mod == "EM" || mod == "FE" || mod == "EF" || 
		   mod == "IE" || mod == "EI" || mod == "EE") {
		# min = E
		PrX1 <- MFIE[4]
		PrX2 <- MFIE[1]
		PrX3 <- MFIE[2]
		PrX4 <- MFIE[3]
	}

	# Probabilities that the conclusion is validated by the n-th test
	# PX1 = 1st test, PX2 = 2nd test, ...
	# i.e., PX1 + PX2 + PX3 + PX4 + residual = 1
	PX1 <- PrX1					# pass
	PX2 <- (1 - PrX1)*PrX2				# fail -> pass
	PX3 <- (1 - PrX1)*(1 - PrX2)*PrX3		# fail*2 -> pass
	PX4 <- (1 - PrX1)*(1 - PrX2)*(1 - PrX3)*PrX4	# fail*3 -> pass

	if (mod == "MM") {
		# min = M
		# M <- PX1; F <- PX2; I <- PX3; E <- PX4
		MFIE <- c(PX1, PX2, PX3, PX4)
	} else if (mod == "MF" || mod == "FM" || mod == "FF") {
		# min = F
		# F <- PX1; I <- PX2; E <- PX3; M <- PX4
		MFIE <- c(PX4, PX1, PX2, PX3)
	} else if (mod == "MI" || mod == "IM" ||
		   mod == "FI" || mod == "IF" || mod == "II") {
		# min = I
		# I <- PX1; E <- PX2; M <- PX3; F <- PX4
		MFIE <- c(PX3, PX4, PX1, PX2)
	} else if (mod == "ME" || mod == "EM" || mod == "FE" || mod == "EF" || 
		   mod == "IE" || mod == "EI" || mod == "EE") {
		# min = E
		# E <- PX1; M <- PX2; F <- PX3; I <- PX4
		MFIE <- c(PX2, PX3, PX4, PX1)
	}

	# The max-heuristic:
	# the most informative premise determines the confidence in
	# the conclusion
	if (mod == "MM" || 
	    mod == "MF" || mod == "FM" ||
	    mod == "MI" || mod == "IM" ||
	    mod == "ME" || mod == "EM") {
		# max = M
		MFIE <- MFIE*um
	} else if (mod == "FF" ||
		 mod == "FI" || mod == "IF" ||
		 mod == "FE" || mod == "EF") {
		# max = F
		MFIE <- MFIE*uf
	} else if (mod == "II" ||
		 mod == "IE" || mod == "EI") {
		# max = I
		MFIE <- MFIE*ui
	} else if (mod == "EE") {
		# max = E
		MFIE <- MFIE*ue
	}

	return(c(MFIE, 1 - sum(MFIE)))
}


#
# To return theoretical probabilities that each syllogistic statement is
# consistent with a Sample Mental Model
#
Syl.simple <- function(mod, fig, x, c, n, m=0.5, f=0.5)
{
	#################################
	# STEP 1 OF THE MODEL
	# To construct a PPM

	prob <- c()
	if (mod == "AA") {
		if (fig == 1) {
			prob <- Syl.ppm.AA1(x, c)	# 1
		}
		else if (fig == 2) {
			prob <- Syl.ppm.AA2(x, c)	# 2
		}
		else if (fig == 3) {
			prob <- Syl.ppm.AA3(x, c)	# 3
		}
		else {	# fig == 4
			prob <- Syl.ppm.AA4(x, c)	# 4
		}
	} else if (mod == "AI" || mod == "AO") {
		if (fig == 1 || fig == 3) {
			prob <- Syl.ppm.AI13(x, c)	# 5a
		}
		else {
			prob <- Syl.ppm.AI24(x, c)	# 6a
		}
	} else if (mod == "IA" || mod == "OA") {
		if (fig == 1 || fig == 2) {
			prob <- Syl.ppm.IA12(x, c)	# 7a
		}
		else {
			prob <- Syl.ppm.IA34(x, c)	# 8a
		}
	} else if (mod == "AE") {
		if (fig == 1 || fig == 3) {
			prob <- Syl.ppm.AE13(x, c)	# 9
		}
		else {
			prob <- Syl.ppm.AE24(x, c)	# 10
		}
	} else if (mod == "EA") {
		if (fig == 1 || fig == 2) {
			prob <- Syl.ppm.EA12(x, c)	# 11
		}
		else {
			prob <- Syl.ppm.EA34(x, c)	# 12
		}
	} else if (mod == "IE" || mod == "OE") {
		prob <- Syl.ppm.IE(x, c)		# 13a
	} else if (mod == "EI" || mod == "EO") {
		prob <- Syl.ppm.EI(x, c)		# 14a
	} else if (mod == "EE") {
		prob <- Syl.ppm.EE(x, c)		# 15
	} else if (mod == "II" || mod == "IO" || mod == "OI" || mod == "OO") {
		prob <- Syl.ppm.II(x, c)		# 16a

	} else if (mod == "AM") {
		if (fig == 1 || fig == 3) {
			prob <- Syl.ppm.AM13(x, c, m)	# 5b
		}
		else {
			prob <- Syl.ppm.AM24(x, c, m)	# 6b
		}
	} else if (mod == "AF") {
		if (fig == 1 || fig == 3) {
			prob <- Syl.ppm.AF13(x, c, f)	# 5c
		}
		else {
			prob <- Syl.ppm.AF24(x, c, f)	# 6c
		}
	} else if (mod == "MA") {
		if (fig == 1 || fig == 2) {
			prob <- Syl.ppm.MA12(x, c, m)	# 7b
		}
		else {
			prob <- Syl.ppm.MA34(x, c, m)	# 8b
		}
	} else if (mod == "FA") {
		if (fig == 1 || fig == 2) {
			prob <- Syl.ppm.FA12(x, c, f)	# 7c
		}
		else {
			prob <- Syl.ppm.FA34(x, c, f)	# 8c
		}
	} else if (mod == "ME") {
		prob <- Syl.ppm.ME(x, c, m)		# 13b
	} else if (mod == "FE") {
		prob <- Syl.ppm.FE(x, c, f)		# 13c
	} else if (mod == "EM") {
		prob <- Syl.ppm.EM(x, c, m)		# 14b
	} else if (mod == "EF") {
		prob <- Syl.ppm.EF(x, c, f)		# 14c
	} else if (mod == "MM") {
		prob <- Syl.ppm.MM(x, c, m)		# 16b
	} else if (mod == "FF") {
		prob <- Syl.ppm.FF(x, c, f)		# 16c
	} else if (mod == "MF") {
		prob <- Syl.ppm.MF(x, c, m, f)		# 16d
	} else if (mod == "FM") {
		prob <- Syl.ppm.FM(x, c, m, f)		# 16e
	} else if (mod == "MI" || mod == "MO") {
		prob <- Syl.ppm.MI(x, c, m)		# 16f
	} else if (mod == "IM" || mod == "OM") {
		prob <- Syl.ppm.IM(x, c, m)		# 16g
	} else if (mod == "FI" || mod == "FO") {
		prob <- Syl.ppm.FI(x, c, f)		# 16h
	} else {# (mod == "IF" || mod == "OF")
		prob <- Syl.ppm.IF(x, c, f)		# 16i
	}

	#################################
	# STEP 2 OF THE MODEL
	# To construct a SMM
	# (calculating theoretical probabilities)

	# The probability that "All S are P" is consistent
	# No sample satisfies (S & not-P)
	A <- (1 - (prob[2] + prob[4]))^n

	# The probability that "Most S are P" is consistent
	# At least one sample satisfies (S & P)
	M <- 1 - (1 - (prob[1] + prob[3]))^n

	# The probability that "Few S are P" is consistent
	# At least one sample satisfies (not-S or not-P)
	F <- 1 - (prob[1] + prob[3])^n

	# The probability that "Some S are P" is consistent
	# At least one sample satisfies (S & P)
	I <- 1 - (1 - (prob[1] + prob[3]))^n

	# The probability that "No S are P" is consistent
	# No sample satisfies (S & P)
	E <- (1 - (prob[1] + prob[3]))^n

	# The probability that "Some S are not P" is consistent
	# At least one sample satisfies (S & not-P)
	O <- 1 - (1 - (prob[2] + prob[4]))^n

	return(c(A, M, F, I, E, O))
}

Allprobs <- function(x, Ps, Pp, P1, P2, P3, P5)
{
	P4 <- Ps - P1 - P2 - P3
	P6 <- x - P1 - P2 - P5
	P7 <- Pp - P1 - P3 - P5
	P8 <- 1 - P1 - P2 - P3 - P4 - P5 - P6 - P7

	pr <- c(P1, P2, P3, P4, P5, P6, P7, P8)
	pr[pr < 0] <- 0

	return(pr)
}

# 1
Syl.ppm.AA1 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- c*x
	Pp <- y

	P1 <- c*x
	P2 <- 0
	P3 <- 0
	P5 <- (1 - c)*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 2
Syl.ppm.AA2 <- function(x, c)
{
	Ps <- c*x
	Pp <- c*x

	P1 <- c*c*x
	P2 <- c*(1 - c)*x
	P3 <- 0
	P5 <- c*(1 - c)*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 3
Syl.ppm.AA3 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- y
	Pp <- y

	P1 <- x
	P2 <- 0
	P3 <- (y - x)*(y - x)/(1 - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 4
Syl.ppm.AA4 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- y
	Pp <- c*x

	P1 <- c*x
	P2 <- (1 - c)*x
	P3 <- 0
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 5a
Syl.ppm.AI13 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- y

	P1 <- x*x
	P2 <- 0
	P3 <- x*(y - x)
	P5 <- x*(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 5b
Syl.ppm.AM13 <- function(x, c, m=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- y

	P1 <- x*x + m*x*(1 - x)
	P2 <- 0
	P3 <- x*(y - x)*(1 - m)
	P5 <- x*(1 - x)*(1 - m)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 5c
Syl.ppm.AF13 <- function(x, c, f=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- y

	P1 <- x*x*(1 - f)
	P2 <- 0
	P3 <- x*(y - x)
	P5 <- x*(1 - (1 - f)*x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 6a
Syl.ppm.AI24 <- function(x, c)
{
	Ps <- x
	Pp <- c*x

	P1 <- c*x*x
	P2 <- (1 - c)*x*x
	P3 <- 0
	P5 <- c*x*(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 6b
Syl.ppm.AM24 <- function(x, c, m=0.5)
{
	Ps <- x
	Pp <- c*x

	P1 <- c*x*(x + m*(1 - x))
	P2 <- (1 - c)*x*(x + m*(1 - x))
	P3 <- 0
	P5 <- c*(1 - m)*x*(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 6c
Syl.ppm.AF24 <- function(x, c, f=0.5)
{
	Ps <- x
	Pp <- c*x

	P1 <- c*(1 - f)*x*x
	P2 <- (1 - c)*(1 - f)*x*x
	P3 <- 0
	P5 <- c*x*(1 - (1 - f)*x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 7a
Syl.ppm.IA12 <- function(x, c)
{
	Ps <- c*x
	Pp <- x

	P1 <- c*x*x
	P2 <- c*x*(1 - x)
	P3 <- 0
	P5 <- (1 - c)*x*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 7b
Syl.ppm.MA12 <- function(x, c, m=0.5)
{
	Ps <- c*x
	Pp <- x

	P1 <- c*x*(x + m*(1 - x))
	P2 <- c*(1 - m)*x*(1 - x)
	P3 <- 0
	P5 <- (1 - c)*x*(x + m*(1 - x))

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 7c
Syl.ppm.FA12 <- function(x, c, f=0.5)
{
	Ps <- c*x
	Pp <- x

	P1 <- c*(1 - f)*x*x
	P2 <- c*x*(1 - (1 - f)*x)
	P3 <- 0
	P5 <- (1 - c)*(1 - f)*x*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 8a
Syl.ppm.IA34 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- y
	Pp <- x

	P1 <- x*x
	P2 <- x*(1 - x)
	P3 <- x*(y - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 8b
Syl.ppm.MA34 <- function(x, c, m=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- y
	Pp <- x

	P1 <- x*x + m*x*(1 - x)
	P2 <- (1 - m)*x*(1 - x)
	P3 <- (1 - m)*x*(y - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 8c
Syl.ppm.FA34 <- function(x, c, f=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- y
	Pp <- x

	P1 <- (1 - f)*x*x
	P2 <- x*(1 - (1 - f)*x)
	P3 <- x*(y - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 9
Syl.ppm.AE13 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- y

	P1 <- 0
	P2 <- 0
	P3 <- x*(y - x)
	P5 <- x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 10
Syl.ppm.AE24 <- function(x, c)
{
	Ps <- x
	Pp <- c*x

	P1 <- 0
	P2 <- 0
	P3 <- 0
	P5 <- c*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 11
Syl.ppm.EA12 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- c*x
	Pp <- x

	P1 <- 0
	P2 <- c*x
	P3 <- 0
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 12
Syl.ppm.EA34 <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- y
	Pp <- x

	P1 <- 0
	P2 <- x
	P3 <- x*(y - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 13a
Syl.ppm.IE <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- 0
	P3 <- x*x
	P5 <- x*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 13b
Syl.ppm.ME <- function(x, c, m=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- 0
	P3 <- (1 - m)*x*x
	P5 <- x*(x + m*(1 - x))

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 13c
Syl.ppm.FE <- function(x, c, f=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- 0
	P3 <- x*x*(1 - (1 - f)*x)/(1 - x)
	P5 <- (1 - f)*x*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 14a
Syl.ppm.EI <- function(x, c)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- x*x
	P3 <- x*x
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 14b
Syl.ppm.EM <- function(x, c, m=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- x*(x + m*(1 - x))
	P3 <- (1 - m)*x*x
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 14c
Syl.ppm.EF <- function(x, c, f=0.5)
{
	y <- 1 - c*(1 - x)
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- (1 - f)*x*x
	P3 <- x*x*(1 - (1 - f)*x)/(1 - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 15
Syl.ppm.EE <- function(x, c)
{
	Ps <- x
	Pp <- x

	P1 <- 0
	P2 <- 0
	P3 <- x*x/(1 - x)
	P5 <- 0

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16a
Syl.ppm.II <- function(x, c)
{
	Ps <- x
	Pp <- x

	P1 <- x*x*x
	P2 <- x*x*(1 - x)
	P3 <- x*x*(1 - x)
	P5 <- x*x*(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16b
Syl.ppm.MM <- function(x, c, m=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- x*(x + m*(1 - x))^2
	P2 <- (1 - m)*x*(1 - x)*(x + m*(1 - x))
#	P3 <- x*x - x*(x + m*(1 - x))^2
	P3 <- (1 - m)^2*x*x*(1 - x)
	P5 <- (1 - m)*x*(1 - x)*(x + m*(1 - x))

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16c
Syl.ppm.FF <- function(x, c, f=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- (1 - f)^2*x*x*x
	P2 <- f*(1 - f)*x*x
#	P3 <- x*x*(1 - (1 - f)^2*x)
	P3 <- x*x*(1 - (1 - f)*x)^2/(1 - x)
	P5 <- f*(1 - f)*x*x

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16d
Syl.ppm.MF <- function(x, c, m=0.5, f=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- (1 - f)*x^2*(x + m*(1 - x))
	P2 <- (1 - m)*(1 - f)*x*x*(1 - x)
#	P3 <- x*x - (1 - f)*x*x*(x + m*(1 - x))
	P3 <- (1 - m)*(1 - (1 - f)*x)*x*x
	P5 <- x*(1 - (1 - f)*x)*(x + m*(1 - x))

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16e
Syl.ppm.FM <- function(x, c, m=0.5, f=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- (1 - f)*x*x*(x + m*(1 - x))
	P2 <- x*(1 - (1 - f)*x)*(x + m*(1 - x))
#	P3 <- x*x - (1 - f)*x*x*(x + m*(1 - x))
	P3 <- (1 - m)*(1 - (1 - f)*x)*x*x
	P5 <- (1 - m)*(1 - f)*x*x**(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16f
Syl.ppm.MI <- function(x, c, m=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- x*x*(x + m*(1 - x))
	P2 <- (1 - m)*x*x*(1 - x)
	P3 <- (1 - m)*x*x*(1 - x)
	P5 <- x*(1 - x)*(x + m*(1 - x))

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16g
Syl.ppm.IM <- function(x, c, m=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- (x + m*(1 - x))*x*x
	P2 <- x*(1 - x)*(x + m*(1 - x))
	P3 <- (1 - m)*x*x*(1 - x)
	P5 <- (1 - m)*x*x*(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16h
Syl.ppm.FI <- function(x, c, f=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- (1 - f)*x*x*x
	P2 <- x*x*(1 - (1 - f)*x)
	P3 <- x*x*(1 - (1 - f)*x)
	P5 <- (1 - f)*x*x*(1 - x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}

# 16i
Syl.ppm.IF <- function(x, c, f=0.5)
{
	Ps <- x
	Pp <- x

	P1 <- (1 - f)*x*x*x
	P2 <- (1 - f)*x*x*(1 - x)
	P3 <- x*x*(1 - (1 - f)*x)
	P5 <- x*x*(1 - (1 - f)*x)

	return(Allprobs(x, Ps, Pp, P1, P2, P3, P5))
}


###############################################################################
#
# Evaluating models
#
# To calculate RMSD for each syllogism 
# between (arcsine transformed [asin(sqrt(p))]) data and model
#
RMSD.syl1 <- function(data1, model1)
{
#	dat <- asin(sqrt(data1))
#	mdl <- asin(sqrt(model1))
#	return((sin(sqrt(sum((mdl - dat)^2)/length(mdl))))^2)
	if (length(model1) == 0) return(0)
	return(sqrt(sum((model1 - data1)^2)/length(model1)))
#	return(sqrt(sum((model1[1:4] - data1[1:4])^2)/4))
}

#
# To calculate Pearson's Correlation Coefficient for each syllogism 
#
CORR.syl1 <- function(data1, model1)
{
#	dat <- asin(sqrt(data1))
#	mdl <- asin(sqrt(model1))
#	return(cor(mdl, dat))
	return(cor(model1, data1))
#	return(cor(model1[1:4], data1[1:4]))
}

#
# To return a vector of RMSDs for all sylogisms 
#
RMSD.syl <- function(data, model, mod=MD)
{
	rmsd <- c()
	lbl <- c()
	for (md in mod) for (fg in c(1:4)) {
		tp <- paste(md, fg, sep="")
		if (is.na(data[tp,1])) next
		rmsd <- c(rmsd, RMSD.syl1(as.matrix(data)[tp,],
					  as.matrix(model)[tp,]))
		lbl <- c(lbl, tp)
	}
	names(rmsd) <- lbl
	return(rmsd)
}


#
# To return a vector of Correlation Coefficients for all sylogisms 
#
CORR.syl <- function(data, model, mod=MD)
{
	corr <- c()
	lbl <- c()
	for (md in mod) for (fg in c(1:4)) {
		tp <- paste(md, fg, sep="")
		if (is.na(data[tp,1])) next
		corr <- c(corr, CORR.syl1(as.matrix(data)[tp,],
					  as.matrix(model)[tp,]))
		lbl <- c(lbl, tp)
	}
	names(corr) <- lbl
	return(corr)
}


#
# for correlation coefficients
#
fisher <- function(r)
{
	r[r > 0.999999999999999] <- 0.9999999999999999
	r[r < -0.999999999999999] <- -0.9999999999999999
	return(0.5*log((1 + r)/(1 - r)))
}


fisher_inv <- function(z)
{
	return((exp(2*z) - 1)/(exp(2*z) + 1))
}


###############################################################################
#
# To return a model (a list of predicted selection rates) 
# besed on parameters "par"
#
PRT.mdl <- function(par, n=7)
{
	rsp <- c(); rn <- c()
	for (md in MD) for (fg in c(1:4)) {
		rn <- c(rn, paste(md, fg, sep=""))
		rsp <- c(rsp, 
			 Syl(md, fg, 
			     par[1], par[2], n, par[3], par[4],
			     par[5], par[6]))
	}
	model <- matrix(rsp, nrow=64, ncol=5, byrow=T)
	rownames(model) <- rn
	colnames(model) <- c("A","I","E","O","N")

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
PRT.par <- function(data, n=7)
{
	#
	# Badness-of-Fit based on 
	# RMSD or Pearson's Correlation Coefficient between data and model
	#
	BOF <- function(x)
	{
		model <- PRT.mdl(x)
		rmsd <- RMSD.syl(data, model)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.40, 0.90,  0.80, 0.60, 0.30, 0.30), BOF,
		method="L-BFGS-B",
		lower=c(0.001, 0.001,  0.001, 0.001, 0.001, 0.001),
		upper=c(0.499, 0.999,  0.999, 0.999, 0.999, 0.999))
	names(rslt$par)=c("x", "c", "ua", "ui", "ue", "uo")

	return(rslt$par)
}


#
# To return a model (a list of predicted selection rates) 
# that best fit to the given data
#
PRT <- function(data, n=7)
{
	param <- PRT.par(data, n)
	model <- PRT.mdl(param, n)

	return(list(p=param, m=model))
}


###############################################################################
#
# To return a model (a list of predicted selection rates) 
# besed on parameters "par"
#
PRT.mdl.AMFO <- function(par, n=7)
{
	rsp <- c(); rn <- c()
	for (md in MD.AMFO) for (fg in c(1:4)) {
		rn <- c(rn, paste(md, fg, sep=""))
		rsp <- c(rsp, 
			 Syl.AMFO(md, fg, 
				  par[1], par[2], n,
				  par[3], par[4], par[5], par[6],
				  par[7], par[8]))
	}
	model <- matrix(rsp, nrow=64, ncol=5, byrow=T)
	rownames(model) <- rn
	colnames(model) <- c("A","M","F","O","N")

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
PRT.par.AMFO <- function(data, n=7)
{
	#
	# Badness-of-Fit based on 
	# RMSD or Pearson's Correlation Coefficient
	# between arcsine transformed data [asin(sqrt(p))] and model
	#
	BOF <- function(x)
	{
		model <- PRT.mdl.AMFO(x)
		rmsd <- RMSD.syl(data, model, MD.AMFO)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model, MD.AMFO)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.40, 0.90,  0.80, 0.60, 0.30, 0.30,  0.50, 0.50), BOF,
		method="L-BFGS-B",
		lower=c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001),
		upper=c(0.499, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999))

	return(rslt$par)
}


#
# To return a model (a list of predicted selection rates) 
# that best fit to the given data
#
PRT.AMFO <- function(data, n=7)
{
	param <- PRT.par.AMFO(data, n)
	model <- PRT.mdl.AMFO(param, n)

	return(list(p=param, m=model))
}


###############################################################################
#
# To return a model (a list of predicted selection rates) 
# besed on parameters "par"
#
PRT.mdl.MFIE <- function(par, n=7)
{
	rsp <- c(); rn <- c()
	for (md in MD.MFIE) for (fg in c(1:4)) {
		rn <- c(rn, paste(md, fg, sep=""))
		rsp <- c(rsp, 
			 Syl.MFIE(md, fg, 
				  par[1], par[2], n,
				  par[3], par[4], par[5], par[6],
				  par[7], par[8]))
	}
	model <- matrix(rsp, nrow=64, ncol=5, byrow=T)
	rownames(model) <- rn
	colnames(model) <- c("M","F","I","E","N")

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
PRT.par.MFIE <- function(data, n=7)
{
	#
	# Badness-of-Fit based on 
	# RMSD or Pearson's Correlation Coefficient
	# between arcsine transformed data [asin(sqrt(p))] and model
	#
	BOF <- function(x)
	{
		model <- PRT.mdl.MFIE(x)
		rmsd <- RMSD.syl(data, model, MD.MFIE)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model, MD.MFIE)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.40, 0.90,  0.80, 0.60, 0.30, 0.30,  0.50, 0.50), BOF,
		method="L-BFGS-B",
		lower=c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001),
		upper=c(0.499, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999, 0.999))
	return(rslt$par)
}


#
# To return a model (a list of predicted selection rates) 
# that best fit to the given data
#
PRT.MFIE <- function(data, n=7)
{
	param <- PRT.par.MFIE(data, n)
	model <- PRT.mdl.MFIE(param, n)

	return(list(p=param, m=model))
}


###############################################################################
#
# Probability Heuristic Model (Chater & Oaksford, 1999)
#
###############################################################################
#
# To return a model (a list of predicted selection rates) 
# besed on parameters "par"
#
PHM.mdl <- function(par)
{
	pa <- par[1]
	pi <- par[2]
	pe <- par[3]
	po <- par[4]
	pent <- par[5]
	perr <- par[6]

	A <- c(rep(pa, 4), rep(perr, 60))
	I <- c(rep(pent, 4), rep(pa, 8), rep(perr, 8), rep(pent, 8), 
	       rep(pi, 4), rep(perr, 8), rep(pent, 8), rep(perr, 4),
	       rep(pent, 12))
	E <- c(rep(perr, 12), rep(pa, 8), rep(perr, 12), rep(pi, 8),
	       rep(perr, 8), rep(pe, 4), rep(perr, 12))
	O <- c(rep(perr, 4), rep(pent, 16), rep(pa, 8), rep(pent, 12),
	       rep(pi, 8), rep(pent, 4), rep(pe, 8), rep(po, 4))
	N <- 1 - (A + I + E + O)

	model <- cbind(A, I, E, O, N)

	rn <- c()
	for (md in MD) for (fg in c(1:4)) {
		tp <- paste(md, fg, sep="")
		rn <- c(rn, tp)
	}
	rownames(model) <- rn
	colnames(model) <- c("A","I","E","O","N")

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
PHM.par <- function(data)
{
	BOF <- function(x)
	{
		model <- PHM.mdl(x)
		rmsd <- RMSD.syl(data, model)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.7, 0.3, 0.2, 0.2, 0.1, 0.1), BOF,
		method="L-BFGS-B",
#		lower=c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
#		upper=c(0.99, 0.99,  0.99, 0.99, 0.99, 0.99))
		lower=c(0.001, 0.001,  0.001, 0.001, 0.001, 0.001),
		upper=c(0.999, 0.999,  0.999, 0.999, 0.999, 0.999))
	names(rslt$par)=c("pa", "pi", "pe", "po", "pent", "perr")

	return(rslt$par)
}


PHM <- function(data)
{
	param <- PHM.par(data)
	model <- PHM.mdl(param)

	return(list(p=param, m=model))
}


###############################################################################
#
PHM.mdl.AMFO <- function(par)
{
	pmin <- par[1]
	pet1 <- par[2]
	pet2 <- par[3]
	perr <- par[4]

	A <- c(rep(pmin, 4), rep(perr, 60))
	M <- c(rep(perr, 4), rep(pmin, 8), rep(perr, 8), rep(pet2, 8), 
	       rep(pmin, 4), rep(perr, 8), rep(pet2, 8), rep(perr, 4),
	       rep(pet2, 12))
	F <- c(rep(perr, 12), rep(pmin, 8), rep(pet2, 8), rep(perr, 4),
	       rep(pmin, 8), rep(pet2, 8), rep(pmin, 4), rep(perr, 12))
	O <- c(rep(perr, 4), rep(pet1, 16), rep(pmin, 8), rep(pet1, 12),
	       rep(pmin, 8), rep(pet1, 4), rep(pmin, 12))
	N <- 1 - (A + M + F + O)

	model <- cbind(A, M, F, O, N)

	rn <- c()
	for (md in MD.AMFO) for (fg in c(1:4)) {
		tp <- paste(md, fg, sep="")
		rn <- c(rn, tp)
	}
	rownames(model) <- rn
	colnames(model) <- c("A","M","F","O","N")

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
PHM.par.AMFO <- function(data)
{
	BOF <- function(x)
	{
		model <- PHM.mdl.AMFO(x)
		rmsd <- RMSD.syl(data, model, MD.AMFO)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model, MD.AMFO)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.70, 0.40, 0.20, 0.05), BOF,
		method="L-BFGS-B",
		lower=c(0.01, 0.01, 0.01, 0.01),
		upper=c(0.99, 0.99, 0.99, 0.99))

	return(rslt$par)
}


PHM.AMFO <- function(data)
{
	param <- PHM.par.AMFO(data)
	model <- PHM.mdl.AMFO(param)

	return(list(p=param, m=model))
}


###############################################################################
#
PHM.mdl.MFIE <- function(par)
{
	pmin <- par[1]
	pet1 <- par[2]
	pet2 <- par[3]
	perr <- par[4]

	M <- c(rep(pmin, 4), rep(perr, 8), rep(pet2, 8), rep(perr, 12), 
	       rep(pet2, 8), rep(perr, 8), rep(pet2, 4), rep(perr, 12))
	F <- c(rep(perr, 4), rep(pmin, 8), rep(pet2, 8), rep(perr, 8), 
	       rep(pmin, 4), rep(pet2, 8), rep(perr, 8), rep(pet2, 4),
	       rep(perr, 12))
	I <- c(rep(pet1, 12), rep(pmin, 8), rep(perr, 8), rep(pet1, 4),
	       rep(pmin, 8), rep(perr, 8), rep(pmin, 4), rep(perr, 12))
	E <- c(rep(perr, 20), rep(pmin, 8), rep(perr, 12), rep(pmin, 8),
	       rep(perr, 4), rep(pmin, 12))
	N <- 1 - (M + F + I + E)

	model <- cbind(M, F, I, E, N)

	rn <- c()
	for (md in MD.MFIE) for (fg in c(1:4)) {
		tp <- paste(md, fg, sep="")
		rn <- c(rn, tp)
	}
	rownames(model) <- rn
	colnames(model) <- c("M","F","I","E","N")

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
PHM.par.MFIE <- function(data)
{
	BOF <- function(x)
	{
		model <- PHM.mdl.MFIE(x)
		rmsd <- RMSD.syl(data, model, MD.MFIE)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model, MD.MFIE)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.40, 0.20, 0.10, 0.10), BOF,
		method="L-BFGS-B",
		lower=c(0.01, 0.01, 0.01, 0.01),
		upper=c(0.99, 0.99, 0.99, 0.99))

	return(rslt$par)
}


PHM.MFIE <- function(data)
{
	param <- PHM.par.MFIE(data)
	model <- PHM.mdl.MFIE(param)

	return(list(p=param, m=model))
}


###############################################################################
#
# Probabilistic Model based on Mental Model Theory
#
###############################################################################
#
# To return a model (a list of predicted selection rates) 
# besed on parameters "par"
#
pMM.mdl <- function(par)
{
	# pn: 第nモデルまで構成する確率

	p1m1 <- par[1]
	p2m1 <- par[2]
	p2m2 <- par[3]
	p3m1 <- par[4]
	p3m2 <- par[5]
	p3m3 <- par[6]

	model <- matrix(nrow=64, ncol=5)
	rn <- c()
	for (md in MD) for (fg in c(1:4)) {
		tp <- paste(md, fg, sep="")
		rn <- c(rn, tp)
	}
	rownames(model) <- rn
	colnames(model) <- c("A","I","E","O","N")

	# 1 Model Syllogisms
	sr1 <- p1m1
	sr2 <- (1 - p1m1)/4

	for (tp in c("AA1", "AA3", "AA4")) {
		model[tp,] <- c(sr1, sr2, sr2, sr2, sr2)
	}
	for (tp in c("AE2", "AE4", "EA1", "EA2")) {
		model[tp,] <- c(sr2, sr2, sr1, sr2, sr2)
	}
	for (tp in c("AI1", "AI3", "IA3", "IA4")) {
		model[tp,] <- c(sr2, sr1, sr2, sr2, sr2)
	}

	# 2 Model Syllogisms
	sr1 <- p2m2
	sr2 <- p2m1
	sr3 <- (1 - p2m1 - p2m2)/3

	for (tp in c("AA2")) {
		model[tp,] <- c(sr2, sr3, sr3, sr3, sr1)
	}
	for (tp in c("EE1", "EE2", "EE3", "EE4",
		     "EO1", "EO2", "EO3", "EO4", "OE1", "OE2", "OE3", "OE4")) {
		model[tp,] <- c(sr3, sr3, sr2, sr3, sr1)
	}
	for (tp in c("AI2", "AI4", "IA1", "IA2", "II1", "II2", "II3", "II4")) {
		model[tp,] <- c(sr3, sr2, sr3, sr3, sr1)
	}
	for (tp in c("AO3", "OA2",
		     "AO2", "OA3",
		     "AO1", "AO4",
		     "OA1", "OA4", 
		     "IO1", "IO2", "IO3", "IO4", 
		     "OI1", "OI2", "OI3", "OI4", 
     		     "OO1", "OO2", "OO3", "OO4")) {
		model[tp,] <- c(sr3, sr3, sr3, sr2, sr1)
	}

	# 3 Model Syllogism
	sr1 <- p3m3
	sr2 <- p3m2
	sr3 <- p3m1
	sr4 <- (1 - p3m1 - p3m2 - p3m3)/2

	for (tp in c("AE1", "AE3",
		     "IE1", "IE2", "IE3", "IE4",
		     "EA3", "EA4",
		     "EI1", "EI2", "EI3", "EI4")) {
		model[tp,] <- c(sr4, sr4, sr3, sr1, sr2)
	}

	return(model)
}


#
# To return an optimal model (parameters) for the data
#
pMM.par <- function(data)
{
	BOF <- function(x)
	{
		model <- pMM.mdl(x)
		rmsd <- RMSD.syl(data, model)
		return(mean(rmsd))
		#corr <- -CORR.syl(data, model)
		#return(fisher_inv(mean(fisher(corr))))
	}

	rslt <-
	  optim(c(0.80, 0.30, 0.50, 0.30, 0.20, 0.30), BOF,
		method="L-BFGS-B",
#		lower=c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
#		upper=c(0.99, 0.99, 0.99, 0.99, 0.99, 0.99))
		lower=c(0.001, 0.001,  0.001, 0.001, 0.001, 0.001),
		upper=c(0.999, 0.999,  0.999, 0.999, 0.999, 0.999))
	names(rslt$par)=c("p1m1", "p2m1", "p2m2", "p3m1", "p3m2", "p3m3")

	return(rslt$par)
}


pMM <- function(data)
{
	param <- pMM.par(data)
	model <- pMM.mdl(param)

	return(list(p=param, m=model))
}


