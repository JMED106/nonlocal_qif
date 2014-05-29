/* Some R utilities */

#include <stdlib.h>
#include <stdio.h>

#include "common.h"
#include "nr.h"
#include "nrutil.h"
#include "Rutils.h"


/* === FUNCTION  R ====================
 * Description:  Calls R script and store variables x, y
 *   Variables:  
 * ======================================= */

void  R_calc(t_data d, double *fr, double *voltage) {
  FILE *Rcmd;
  char cmd[30], cmd2[30], cmd3[30], cmd4[30];
  Rcmd = popen("R -e \"source('./volt_avg.R')\" -q -s","r");
  fscanf(Rcmd,"%s %s\n%s %s",cmd, cmd2,cmd3,cmd4);
  pclose(Rcmd);
  *fr =atof(cmd2);
  *voltage = atof(cmd4);
}


/* === FUNCTION  R ====================
 * Description:  Creates R script for R_calc
 *   Variables:  vpeak, aprox x and y (theta)
 * ======================================= */

void R_script(t_data d, double x, double y) {
  FILE *Rscript;
  Rscript = fopen("volt_avg.R","w");
  
  fprintf(Rscript,"setwd('./temp');\n"
	  "volt_dist_files = list.files(pattern= 'volt_dist*');\n"
	  "voltages <- numeric(length(volt_dist_files));\n"
	  "firing_rates <- numeric(length(volt_dist_files));\n"
	  "for (i in volt_dist_files) \{\n"
	  "    data<-read.table(i);\n"
	  "    data2 <- data[abs(data[,1])<=%d,];\n"
	  "    data3 <- hist(data2,breaks=seq(-%d,%d,by=0.1),freq=FALSE);\n"
	  "    x <-data3$breaks[1:length(data3$counts)]\n"
	  "    datax <- data.frame(x,data3$density)\n"
	  "    lorentz <- nls(data3$density ~(1/pi)*a/((x-b)*(x-b)+a*a),data=datax,start=list(a=%lf,b=%lf),algorithm=\"plinear\", nls.control(maxiter = %d , tol = %e, minFactor = 1/%d, printEval = FALSE, warnOnly = FALSE))\n"
	  "    voltages[i] <- coef(lorentz)[2];\n"
	  "    firing_rates[i] <- coef(lorentz)[1];\n"
	  "}\n"
	  "avg_volt <- mean(voltages);\n"
	  "avg_fr <- mean(firing_rates);\n"
	  "print(avg_fr);\n"
	  "print(avg_volt);\n",(int)d.vp,(int)d.vp,(int)d.vp,x*(M_PI),y,50,1e-5,2048);
  fclose(Rscript);
}
