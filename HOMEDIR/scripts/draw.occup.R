draw.occup <- function(analysis = "wt", center = 3109, window = 250, ymax.occup = 1.1, 
			wt.length = 6001, 
			annotation = data.frame(name = "", color = "", left = 0, 
            		right = 0)[0,]){
	
	ylim <- c(-0.1, ymax.occup)

	result <- get(analysis)
	seqlength <- nrow(result)/5
	ori.center <- center
	center <- center + seqlength*2
	insert.length <- seqlength - wt.length

	plot(x = result$pos[(center - window):(center + window)], 
		y = result$nucoccup[(center - window):(center + window)], 
		xlim = c(ori.center-window, ori.center+window), 
		ylim = ylim, type = "n", lwd = 2, col = "blue", xaxt = "n", yaxt = "n", 
		xlab = "Position (bp)", ylab = "Occupancy", main = analysis)
	polygon(x = c(ori.center-window, (ori.center-window):(ori.center+window), ori.center+window), 
		y = c(0, result$nucoccup[(center-window):(center+window)], 0), col = "gray")
	box(lwd = 2)
	axis(side = 1, lwd.ticks = 2, at = seq(ori.center-300, ori.center+300, 100), 
			labels = c(-300, -200, -100, 0, 100, 200, 300))
	axis(side = 2, lwd.ticks = 2, at = seq(0, 1, 0.5), 
			labels = seq(0, 1, 0.5), las = 2)
	
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

        for(i in seq_len(nrow(annotation))){
            xleft <- annotation$left[i]
            xright <- annotation$right[i]
            color <- annotation$color[i]
            rect(xleft = xleft, ybottom = -0.08, xright = xright, ytop = -0.05, 
                col = color, border = color)
        }

	cat("insert.length: ", insert.length, "\n")
        for(i in seq_len(nrow(annotation))){
		cat(annotation[i,1], "\t", annotation[i,2], "\t", 
			annotation[i,3], "\t", annotation[i,4], "\t", "\n")
	}
}
