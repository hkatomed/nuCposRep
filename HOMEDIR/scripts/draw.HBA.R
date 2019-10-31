draw.HBA <- function(analysis = "wt", center = 3109, window = 250, ylim = c(-35, 15), 
			wt.length = 6001, ylabels = c(-20, 0, 20), 
			annotation = data.frame(name = "", color = "", left = 0, 
            		right = 0)[0,], ablines = list(v=0, h=0), 
			col.ablines = list(v="gray", h="gray"), 
			lty.ablines = list(v=1, h=1)){
	result <- get(analysis)
	seqlength <- nrow(result)/5
	ori.center <- center
	center <- center + seqlength*2
	insert.length <- seqlength - wt.length

	plot(x = result$pos[(center - window):(center + window)], 
		y = result$affinity[(center - window):(center + window)], 
		xlim = c(ori.center-window, ori.center+window), 
		ylim = ylim, type = "n", lwd = 2, col = "blue", xaxt = "n", yaxt = "n", 
		xlab = "Position (bp)", ylab = "HBA", main = analysis)
	abline(v = ori.center + ablines$v, lwd = 2, col = col.ablines$v, lty = lty.ablines$v)
	abline(h = ablines$h, lwd = 2, col = col.ablines$h, lty = lty.ablines$h)
	lines(x = result$pos[(center - window):(center + window)], 
		y = result$affinity[(center - window):(center + window)], 
			lwd = 2, col = "blue")
	box(lwd = 2)
	axis(side = 1, lwd.ticks = 2, at = seq(ori.center-300, ori.center+300, 100), 
			labels = c(-300, -200, -100, 0, 100, 200, 300))
	axis(side = 2, lwd.ticks = 2, at = ylabels, labels = ylabels, las = 2)
	
	if(insert.length > 0){
		anno.insertion = data.frame(name = "ins", color = "red", left = ori.center, 
            		right = ori.center + insert.length - 1)
		annotation <- rbind(annotation, anno.insertion)
	}

	if(nrow(annotation) > 0){
		for(i in 1:nrow(annotation)){
			if(annotation$left[i] > ori.center){
				annotation$left[i] <- annotation$left[i] + insert.length
				annotation$right[i] <- annotation$right[i] + insert.length
			}
			if(annotation$left[i] < ori.center-window){
				annotation$left[i] <- ori.center+window
			}
			if(annotation$right[i] > ori.center+window){
				annotation$right[i] <- ori.center+window
			}
		}
	}

	ybottom <- ylim[1] + 1
	ytop <- ylim[1] + 3
        for(i in seq_len(nrow(annotation))){
            xleft <- annotation$left[i]
            xright <- annotation$right[i]
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = ybottom, xright = xright, 
                ytop = ytop, 
                col = color, border = color)
        }

	cat("insert.length: ", insert.length, "\n")
        for(i in seq_len(nrow(annotation))){
		cat(annotation[i,1], "\t", annotation[i,2], "\t", 
			annotation[i,3], "\t", annotation[i,4], "\t", "\n")
	}

}
