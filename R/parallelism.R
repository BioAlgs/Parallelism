#' Interaction test for bivariate smoothing spline ANOVA
#'
#' @import gss
#' @param y the input signal.
#' @param z the spatial or tempral locations (scaled to [0,1]).
#' @param g the categorical variable representing the group information.
#'
#' @return score: the value of the test statisitcs.
#' @return pvalue: P-value of the test.
#'
#' @examples
#' \donttest{
#' library(gss)
#' library(Parallelism)
#' n = 200
#' ##generate locataions of the spatial or temperal points
#' z= seq(0,1,length.out=n)
#' z = c(z,z)
#' ##the group indictor for each mearsuement
#' g = c(rep(0,n),rep(1,n))
#' df= data.frame(g=as.factor(g),z=z)
#' ##generate the measured signal for each group at different locations
#' fun1 = function(x,delta) 2.5*sin(3*pi*x[,2])*(1-x[,2])*ifelse(x#' [,1]==0,1,0) + ( (2.5+delta)*sin(3*pi*x[,2])*(1-x[,2]) )*#' ifelse(x[,1]==0,0,1)
#' y = fun1(df, 1) + rnorm(2*n)
#' 
#' ##begin the test
#' res = parallelism(y,z,as.factor(g))
#'}
#' @export parallelism
#' @import gss

parallelism=function(y=y, z=z, g=g){

    ###Check the input data
    g = as.factor(g)
    groupName = levels(g)
    if(length(groupName)==1){
        stop('There is only one group.')
    }

    site = sort(z[g==groupName[1]])
    for(i in 2:length(groupName)){
        if(!all(site==sort(z[g==groupName[i]]))){
            stop('The location variable in each are not consistent')
        }
    }

    ###fitting model
    cat('begin model fitting...\n')
    df = data.frame(y=y, g=as.factor(g), z=z)
    ps = ssanova0(y ~ g * z ,type=list(z=list("per",c(0,1))), data=df)
    cat('begin testing...\n')

    ###begin testing 
        object = ps
        newdata=df
        # fitted(ps)
        newdata=ps$mf
        include=c(object$terms$labels,object$lab.p)
        # plot(df$z, fitted(ps))
        # predict(ps, newdata, include=c("g:z"))


        mf=ps$mf
        mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
        mf$prec <- mf$maxiter <- NULL
        nnew <- dim(newdata)[1]
        nobs <- length(object$c)
        nnull <- length(object$d)
        labels.p <- object$lab.p
        ## Extract included terms
        term <- object$terms
        philist <- rklist <- NULL
        s <- q <- NULL
        nq <- 0
        for (label in include) {
            if (label=="1") {
                philist <- c(philist,term[[label]]$iphi)
                s <- cbind(s,rep(1,len=nnew))
                next
            }
            if (label%in%labels.p) next
            if (label=="offset") next
            xnew <- newdata[,term[[label]]$vlist]
            x <- object$mf[,term[[label]]$vlist]
            nphi <- term[[label]]$nphi
            nrk <- term[[label]]$nrk
            if (nphi) {
                iphi <- term[[label]]$iphi
                phi <- term[[label]]$phi
                for (i in 1:nphi) {
                    philist <- c(philist,iphi+(i-1))
                    s <- cbind(s,phi$fun(xnew,nu=i,env=phi$env))
                }
            }
            if (nrk) {
                irk <- term[[label]]$irk
                rk <- term[[label]]$rk
                for (i in 1:nrk) {
                    rklist <- c(rklist,irk+(i-1))
                    nq <- nq+1
                    q <- array(c(q,rk$fun(xnew,x,nu=i,env=rk$env,out=TRUE)),c(nnew,nobs,nq))
                }
            }
        }
        qq <- matrix(0,nnew,nobs)
        nq <- 0
        for (i in rklist) {
            nq <- nq + 1
            qq <- qq + 10^object$theta[i]*q[,,nq]
        }


        K2 = 10^object$theta[2]*q[,,2]
        pmean_inter = K2%*%object$c
        M = qq + 10^object$nlambda * diag(rep(1, nnew))
        delta = K2 %*%  (solve(M) - solve(M) %*% s %*% solve(t(s) %*% solve(M) %*% s) %*% t(s) %*% solve(M))

        etahat =delta %*% object$mf[,1]
        tn = sum((delta %*% object$mf[,1])^2)
        Delta = t(delta) %*% delta
        mu = sum(diag(Delta))
        sigma = sqrt(2*sum(diag(Delta%*%Delta)))
        score = abs(tn - mu)/sigma
    pval = 2*pnorm(-abs(score))
    cat('the test complete.\n')
    cat(paste0('p-value=',pval),'\n')
    return(list(score=score, pvalue=pval))
}
