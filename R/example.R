fun1 = function(x,delta) 2.5*sin(3*pi*x[,2])*(1-x[,2])*ifelse(x[,1]==0,1,0) + ( (2.5+delta)*sin(3*pi*x[,2])*(1-x[,2]) )* ifelse(x[,1]==0,0,1)


png("curve1.png",width = 480, height = 360)
curve_df = cbind(fun1(df, 0.5)[1:n], fun1(df, 1)[(n+1):(2*n)])
matplot(curve_df, type="l", pch=1,xlab="time",ylab="signal",xaxt="n")
axis(1, at=seq(0,100,20), label=seq(0,1,length.out=length(seq(0,100,20))))
legend("topright", legend = c("group 0","group 1"), col=1:2, lty=1:2)
dev.off()
