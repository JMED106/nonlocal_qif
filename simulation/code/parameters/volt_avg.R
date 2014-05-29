setwd('./temp');
volt_dist_files = list.files(pattern= 'volt_dist*');
voltages <- numeric(length(volt_dist_files));
firing_rates <- numeric(length(volt_dist_files));

for (i in volt_dist_files) {
    data<-read.table(i);
    data2 <- data[abs(data[,1])<=50,];
    data3 <- hist(data2,breaks=seq(-50,50,by=0.1),freq=FALSE);
    x <-data3$breaks[1:length(data3$counts)]
    datax <- data.frame(x,data3$density)
    lorentz <- nls(data3$density ~(1/pi)*a/((x-b)*(x-b)+a*a),data=datax,start=list(a=1,b=0))
    voltages[i] <- coef(lorentz)[2];
    firing_rates[i] <- coef(lorentz)[1];
}
avg_volt <- mean(voltages);
avg_fr <- mean(firing_rates);
print(avg_fr);
print(avg_volt);