############################# PMM.R ##############################
#
# PMM. R   -    Predictive Maintenance Module, v1.0
#
#               Input Data: Ball-Screw Monitored Time Series: 5th
#
#               Florian Sobieczky - Nov 2022 - SCCH
#
# Usage:
#
#   D <- agg(direc='Path/', I=2:5, feats=features1, )
#
# This aggregates the features in vector 'feats' over the cycles and 
# stores them in D.
#
#   show_D(D)
#
# to see it again. Then do the anomaly-detection, by picking one
# of the features - in this case column 3 is pretty good (7, 9, and
# 15 also!), and say
#
#   D <- detect_anomaly(D, feat=3, perc_thresh=10)
#
# In verbose (=default) mode, this shows the standardized feature 
# and at which it surpasses the threshold-level at least perc_thresh
# percent of the time.
#
# Now, D contains a label D$hi, with which training of
# ML models for health-state classification can be carried out.
#
#######################################################################

library(data.table)
library(tsfeatures)
library(qcc)
library(caret)
library(rpart)

all_features <- c("entropy", "lumpiness", "stability", "max_level_shift", "max_var_shift", "max_kl_shift", "crossing_points", "flat_spots", "hurst", "unitroot_kpss", "unitroot_pp", "stl_features", "acf_features", "pacf_features", "holt_parameters", "hw_parameters", "heterogeneity", "nonlinearity", "arch_stat", "compengine", "embed2_incircle", "ac_9", "firstmin_ac", "firstzero_ac", "trev_num", "motiftwo_entro3", "binarize_mean", "walker_propcross", "localsimple_taures", "sampen_first", "sampenc", "std1st_der", "spreadrandomlocal_meantaul", "histogram_mode", "outlierinclude_mdrmd", "fluctanal_prop_r1")

features1 <- c("entropy", "lumpiness", "stability", "max_level_shift", "max_var_shift")


########################### Reading in Data ##########################################


# Change format of Timestamp

fix_time <- function(D, time_col=which(names(D)=="Timestamp")) {
    n<-dim(D)[1]

    fixed_t <- rep('', n)
    for(i in 1:n) {
        t <- D[i, time_col]
        s <- strsplit(as.character(t), ":")[[1]]
        ft <- paste0(s[1], ':', s[2], ':', s[3], '.', s[4])

        s <- strsplit(ft, " ")[[1]]
        ft <- paste0(s[1], '-', s[2], '-', s[3], ' ', s[4])
        fixed_t[i]<- ft
    }
    D$TS <- fixed_t  

    D
}


# Compose the filename with path.

get_fname <- function(i, direc='../data/Fifth/raw/') {    
    files <- dir(direc)                #  Reading in all filenames
    fname <- paste0(direc, files[i])

    fname
}


# 

read_A <- function(fname, num_parts=2, i_part=1) {
       A <- fread(fname)                 # If data.table not available, use read.csv():
       # A <- read.csv(fname, sep=';', skip=4)    # fread() from lib data.table faster!
       A <- as.data.frame(A)        
       A <- fix_time(A)                  # Reformulation of time + Date.
       for(i in 1:length(names(A))) {
           names(A)[i] <- gsub("\\.+", "_", names(A)[i])
           names(A)[i] <- gsub(" ", "_", names(A)[i])
       }

       # Now divide into parts parts and select part n_part:

       m <- floor(dim(A)[1]/num_parts) 
       A <- A[((i_part-1)*m+1):(i_part*m),]

       # Return finished data.frame:
       A
}


########################### Aggregating Data over cycles ##################

agg <- function(direc='../data/Fifth/raw/', num_parts=2, i_part=1, I=2:5, feats=features1, verbose=TRUE) {
    # cycle_parts: Number of homogeneous parts of single cycle: Here 2 (for- and backward)
    # I: vector of columns which contain raw features to be considered in Anomaly Detection

    n_raw <- length(I)                   # number of raw features (measured sensor values)
    m <- length(feats)                   # number of derived features 
    n_all <- m * n_raw                   # number of all features
    files <- dir(direc)                  # Reading in all filenames
    n <- length(files)                   # number of files = number of lines in D
    D <- matrix(rep(0, n*(n_all+1)), ncol=n_all+1) # Initializes Result matrix
    D <- as.data.frame(D)                # D is a data frame.
    
    for(i in 1:n) {
        if(verbose==TRUE)
            print(paste(i,'from',n,' : m=', m, ' -  n_raw=', n_raw))
        fname <- get_fname(i, direc=direc)             # Name of current file.
        A <- read_A(fname, num_parts=num_parts, i_part=i_part)
                                        # A is dataframe of current file.

        G <- data.frame()
        G[1,1] <- A$TS[1]                    # Take time variable (begin of cycle)
        for(j in 1:length(I)) {
            F <- as.data.frame(tsfeatures(A[,I[j]], features=feats))
            G[1,(1+(j-1)*m+1):(1+j*m)] <- as.numeric(F[1,1:m])
        }

        D[i, ] <- G
    }

    namesD <- rep('', n_all)   # Naming of the columns: "derived_feature(raw_feature)"
    namesD[1] <- 'TS'
    for(j in 1:n_raw)
        for(l in 1:m)
            namesD[(j-1)*m+l+1] <- paste0(feats[l],'(',names(A)[I[j]],')')
    names(D) <- namesD
    
    D
}


################################# Plotting ############################################

show_D <- function(D, b=3, l=2, wait=TRUE) {
    par(mfrow=c(l, b))

    n_all <- dim(D)[2]-1                        # Total number of features (=m*n_raw)
    n_pics <- ceiling(n_all/(b*l))              # Number of plots (with multiple diagrams!)
    n_rest <- b*l-(n_pics*b*l - n_all)          # Num of diagrams on last plot.

    for(k in 1:(n_pics-1)) {
        for(i in ((k-1)*b*l+2):(k*b*l+1)){
            plot(as.POSIXct(D$TS), D[,i], main=names(D)[i], xlab='Time', ylab='', t='l')
            grid()
        }
        if(wait==TRUE)
            readline()
    }
    print(n_rest)
    readline()
    if(n_rest>0)
        for(i in ((n_pics-1)*b*l+2):(n_all+1)) {
            plot(as.POSIXct(D$TS), D[,i], main=names(D)[i], xlab='Time', ylab='', t='l')
            grid()
        }
        
}


############################# Detection using Predict & Compare ########################

detect_anomaly <- function(D, feat=3, fun=mean, perc_thresh=10, lb_size=40, pr_size=20, th_grid=seq(0.3,4,length.out=100), verbose=TRUE) { 

    n <- dim(D)[1]       # Length of considered interval (of vector of aggregated data).
    n_all <- dim(D)[2]   # 

    
    v_m <- mean(as.numeric(D[,feat]))
    v_s <- sd(as.numeric(D[,feat]))
    v <- standardize(as.numeric(D[,feat]))        # This vector to assign labels.

    comp_v <- pa_v <- pr_v <- rep(0, n)         # past(lookback) and prediction vectors.
    
    for(i in (lb_size+1):(n-pr_size)) {
        input_vec <- v[(i-lb_size):(i-1)]
        present_vec <- v[i:(i+pr_size)]

        pa_v[i] <- past_value <- fun(input_vec)
        pr_v[i] <- present_value <- fun(present_vec)

        comp_v[i] <- pr_v[i] - pa_v[i]          # This is the health indicator.
    }

    comp_v[(n-pr_size+1):n] <- comp_v[n-pr_size]
    D$hi_feature <- comp_v

    A <- data.frame(D$TS, v, pa_v, pr_v, comp_v)
    names(A) <- c('TS', 'v', 'pa_v', 'pr_v', 'comp_v')

    perc_above <- rep(100, length(th_grid))
    for(j in 1:length(th_grid)) {
        s <- th_grid[j]
        n_above <- length(which(comp_v > s))
        if(length(n_above)==0)
            n_above <- 0
        perc_above[j] <- 100*n_above/n
    }
    i_ideal <- min(which(perc_above < perc_thresh))
    print(perc_above)
    print(i_ideal)
    print(th_grid[i_ideal])
    final_threshold <- un_standardize(th_grid[i_ideal], v_m, v_s)
    
    D$hi_feature <- comp_v
    D$hi <- ifelse(comp_v>th_grid[i_ideal], 1, 0)

    if(verbose==TRUE) {
        show_A(A, title=paste('HI:', names(D)[feat],'   /  Thresh = ', round(th_grid[i_ideal], 2)) , thresh=th_grid[i_ideal])
        readline()
        show_finished_D(D, feat, perc_threshold[i_ideal], final_thresh=round(final_threshold, 2), title=paste('HI: ', names(D)[feat],'  /  Thresh = ', round(final_threshold, 2)) )
    }

        
    D
}


show_A <- function(A, title='', thresh=2.5) {
    par(mfrow = c(1,1))

    plot(as.POSIXct(A$TS), A$v, main=title, xlab='Time', ylab='', t='l')
    lines(as.POSIXct(A$TS), A$comp_v, col='red')
    abline(h=thresh, col='red', lt=2)
    grid()    
}


show_finished_D <- function(D, feat, perc_thresh=perc_thresh, final_threshold=0, title='') {
    plot(as.POSIXct(D$TS), D[,feat], t='l', main=title, xlab='Time', ylab=names(D)[feat])
    lines(as.POSIXct(D$TS), D$hi, t='l', col='blue', lwd=3)
    abline(h=final_threshold, col='red', lt=2)
    grid()
}

standardize <- function(v) (v-mean(v))/sd(v)

un_standardize <- function(v, m, s) v*s+m


############################## Classify State of Health ##########################


train_mod <- function(D, feat, train_interval=c(1:200), lb_size=40, noise=0.3) {

    v_m <- mean(D[, feat])
    v_s <- sd(D[, feat])
    v <- standardize(D[, feat])
    v_new <- v + rnorm(length(v), 0, noise)
    feat_new <- un_standardize(v_new, v_m, v_s)    

    D0 <- data.frame(D$TS, D[,feat], D$hi)
    D0$hi <- as.factor(D$hi, levels=c(0,1))
    Dnew <- data.frame(D$TS, feat_new, D$hi)

    names(D0) <- names(D) <- c('Time', 'feature', 'hi')

#    mod <- train(hi~., data=D0, method="rf", metric="RMSE")
    mod <- rpart(hi~., data=D0, method="class")
    pred <- predict(mod, newdata=Dnew)
    
    list(mod, pred)
}



############################## Prognose RUL ##################################

