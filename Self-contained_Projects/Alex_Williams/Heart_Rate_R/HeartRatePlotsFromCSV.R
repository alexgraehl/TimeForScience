require("stats") # for supsmu (super smooth)

# This is for reading a CSV file of heart rate data as exported by the "Wahoo" app for the iPhone / iPad
options(stringsAsFactors=FALSE)
options(error=recover)

WAHOO_FITNESS_DIRECTORY <- "~/Dropbox/Apps/WahooFitness"

zips <- list.files(path=WAHOO_FITNESS_DIRECTORY, recursive=T, pattern="*.zip", full.names=T)
for (f in zips) {
	csvName <- sub("zip$","csv", f, perl=T, ignore.case=T)
	if (file.exists(csvName)) {
		print(paste("CSV file <", csvName, "> appears to already have been extracted. Skipping unzip.", sep=''))
	} else {
		unzip(f, exdir=dirname(f));
		#system(paste("bzip2 ", csvName, sep=''))
		#stopifnot(file.exists(csvName))
	}
}

csvs <- list.files(path=WAHOO_FITNESS_DIRECTORY, recursive=T, pattern="*.csv", full.names=T)

for (fff in csvs) {
	#"~/Desktop/2012-06-03_0101_Golfing_WF.csv" #HRM 5-20 after is when slpn/2012-06-02_0205_Golfing_WF.csv"
	
	#system(
	
	FIRST_INTERESTING_HEADER_REGEXP <- "^time,pwr_accdist,pwr_cadence" # regular expression
	MAX_LINES_TO_LOOK_FOR_INTERESTING_HEADER <- 500
	findHeader <- readLines(fff, n=MAX_LINES_TO_LOOK_FOR_INTERESTING_HEADER) #, blank.lines.skip=F)
	header.vec <- grep(FIRST_INTERESTING_HEADER_REGEXP, findHeader, perl=T, ignore.case=T)
	if (length(header.vec) != 1) {
		warning(paste("ERROR: No single header line was found in the file <", fff, ">. Maybe it's not the right type of CSV file. Skipping."))
		next;
	}
	print(paste("Reading the file <", fff, ">...", sep=''))
	
	pdfname <- sub("csv", "pdf", fff, perl=T, ignore.case=T)
	pngname <- sub("csv", "png", fff, perl=T, ignore.case=T)
	
	## For filtering out totally bizarre broken values from the data
	MAX_THEORETICAL_HEARTRATE <- 500
	MIN_THEORETICAL_HEARTRATE <-  25
	# ========================================
	NUM_USELESS_TOP_LINES <- (header.vec[1] - 1) ## num lines in the input file before we get to the table part
	dat.all <- read.csv(fff, skip=NUM_USELESS_TOP_LINES, header=T)
	# ========================================
	
	dd <- dat.all[, c("time","hr_heartrate")] ## only get the fields of interest
	
	totalSeconds <- as.integer(max(dd[, "time"]))
	totalMinutes <- totalSeconds/60
	totalHours   <- totalMinutes/60
	
	hr.unfiltered <- dd[, "hr_heartrate"]
	
	hr <- hr.unfiltered
	hr[ hr > MAX_THEORETICAL_HEARTRATE | hr < MIN_THEORETICAL_HEARTRATE ] <- NA  # remove totally broken values from the data
	hintORIGINAL <- as.integer(round(hr)) #as.integer(round(hr))
	
	filterOutLongRunsOfSameValue <- function(windowSize=30, input.vec) {
		# if heart rate doesn't change by even one beat for this many seconds, then this is probably bad data, so NA it out
		result.vec <- input.vec
		for (i in seq(windowSize, length(input.vec))) {
			if (is.na(input.vec[i])) {
				next;
			} else {
				whichIndices <- seq(i-windowSize+1, i)
				vals <- input.vec[whichIndices]
				theyAreNA <- is.na(vals)
				notNAAndSame <- (vals==input.vec[i]) & (!is.na(vals))
				if (all(theyAreNA | notNAAndSame)) {
					result.vec[i] <- NA
				}
			}
		}
		return(result.vec)
	}
	
	hint <- filterOutLongRunsOfSameValue(windowSize=30, hintORIGINAL)
	
	htab <- cbind(seq_along(hint), hint)
	colnames(htab) <- c("seconds", "bpm")
	
	hfinite <- htab[ is.finite(htab[,"bpm"]), ]
	beats <- hfinite[,"bpm"]
	sec <- hfinite[,"seconds"]
	minutes <- sec/60
	
	numSecondsToIntegrateForSmoothing <- 30
	desiredSpanFraction <- numSecondsToIntegrateForSmoothing/totalSeconds # smooth to exactly one minute worth of data
	
	ss <- stats::supsmu(x=minutes, y=beats, span=desiredSpanFraction) # super smooth
	
	removeGapsFromSupSmu <- function(xGapMustBeThisBigInMinutes, inputSS) {
		checkUs <- seq(from=2, to=length(inputSS$x))
		TOO_BIG_X_GAP <- 1 
		tooBigGap.bool.vec <- (inputSS$x[checkUs] - inputSS$x[checkUs-1]) > xGapMustBeThisBigInMinutes
		inputSS$x[tooBigGap.bool.vec] <- NA
		return(inputSS)
	}
	
	ss <- removeGapsFromSupSmu(xGapMustBeThisBigInMinutes=1.00, ss) # gap of 1 minute or more = this is a big gap, show it as a gap in the data!
	maxValueInSmoothedLine <- max(ss$y, na.rm=T)
	
	JAGGED_LINE_COL <- "#FF000044"
	SMOOTH_LINE_COL <- "#000000FF"
	LOW_HR_COL <- "#9999CC66"
	MED_HR_COL <- "#9999CC33"
	HIGH_HR_COL <- "#FFCC0044"
	HOUR_BAR_COL <- "#00000066"
	
	LOW_HR_AT <- 55     # Where the "low heart rate" region is drawn
	HI_HR_AT  <- 85     # Where the "high heart rate" region is drawn
	MIN_PLOT_YVAL <- 35 # Lowest y-value in the plot
	MAX_PLOT_YVAL <- 120 # Highest y-value in the plot
	
	PIXELS_PER_HOUR_IN_PNG <- 400
	PNG_PLOT_HEIGHT <- 700
	
	## ============== Generate the plot ======================	
	pngWidth <- max(500, totalHours*PIXELS_PER_HOUR_IN_PNG)
	
	pdfWidth <- pngWidth / 100
	pdfHeight <- PNG_PLOT_HEIGHT / 100
	
	#png(pngname, width=pngWidth, height=PNG_PLOT_HEIGHT, pointsize=18, res=72)
	pdf(pdfname, width=pdfWidth, height=pdfHeight, pointsize=12)
	
	par(las=1) # always horizontal axis labels!
	#par(xaxp=c(0, totalMinutes, 100))
	plot(ss, type='n', main=paste("Heart Rate: BPM over", round(totalHours,1), "hours")
		, ylab="Beats per minute"
		, xlab="Minutes"
		, xaxt='n'
		, ylim=c(MIN_PLOT_YVAL, MAX_PLOT_YVAL)
		) # blank plot!
		
	xlabels <- seq(from=0, to=totalMinutes+1, by=10)
	hourstarts <- seq(from=0, to=floor(totalHours)+1, by=1) # in minutes
	everyOtherHourStarts <- seq(from=0, to=floor(totalHours)+1, by=2) # in minutes
	axis(side=1, at=xlabels, labels=xlabels) #, pos=, lty=, col=, las=, tck=, ...)
	
	rect(xleft=(-totalMinutes), xright=totalMinutes*2, ybottom=0, ytop=LOW_HR_AT, col=LOW_HR_COL, border=NA)
	#rect(xleft=(-totalMinutes), xright=totalMinutes*2, ybottom=50, ytop=65, col=MED_HR_COL, border=NA)
	rect(xleft=(-totalMinutes), xright=totalMinutes*2, ybottom=HI_HR_AT, ytop=500, col=HIGH_HR_COL, border=NA)
	
	rect(xleft=everyOtherHourStarts*60, xright=everyOtherHourStarts*60+60, ybottom=0, ytop=MIN_PLOT_YVAL+5, col=HOUR_BAR_COL, border=NA)
	#text(x=hourstarts*60, y=MIN_PLOT_YVAL+2, seq_along(hourstarts), pos=4)

	
	RANDOM_JITTER_FACTOR <- rnorm(length(beats), sd=0.50)
	points(x=minutes, y=beats+RANDOM_JITTER_FACTOR, col=JAGGED_LINE_COL, pch='.', cex=3.5)
	lines(ss, col=SMOOTH_LINE_COL, lwd=2)
	
	dev.off()
	
	## =================== Done making a plot ======================
	
	# ==============================================================
	
	#hfinite.df <- data.frame("seconds"=hfinite[,"seconds"], "bpm"=hfinite[,"bpm"])
	#require("ggplot2")
	#zzz <- ggplot2::ggplot(hfinite.df, aes(seconds, bpm)) ## CANNOT quote the things in "aes" or it breaks
	#zzz + geom_point() + stat_smooth() # makes a plot
	
	# ==============================================================
}
	
