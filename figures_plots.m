
%%%%%%%%%%%%%%%%%%% Plot the degree distributions of all competetive models
%%%%%%%%%%%%%%%%%%% correponding to Twitter Data Set.
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%

clear all;


data=csvread('output_twitter.csv',1);
xDeg=data(:,1);
xAct=data(:,2);
xlom=data(:,3);
xver1=data(:,4);
xver2=data(:,5);
xver3=data(:,6);
xver4=data(:,7);
xPow=data(:,9);
xPar=data(:,10);
xLog=data(:,11);
xPoC=data(:,12);
xExp=data(:,13);

figure

loglog(xDeg,xAct,'.','MarkerSize',7,'MarkerEdgeColor','b')
%semilogy(xDeg,xAct,'.','MarkerSize',7)
%plot(xDeg,xAct,'.','MarkerSize',7)

%txt1 = 'Actual';
%text(xDeg(10),xAct(10),txt1)

hold on
loglog(xDeg,xPar,'linewidth',1,'color',[0, 0.75, 0.75])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xPow,'linewidth',1,'color',[0.4940, 0.1840, 0.5560])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xLog,'linewidth',1,'color',[0.75,0,0.75])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)


hold on
loglog(xDeg,xPoC,'linewidth',1,'color',[0.75,0.75,0])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xExp,'linewidth',1,'color',[0.25,0.25,0.25])
%semilogy(xDeg,xEst,'linewidth',1.5)
%plot(xDeg,xEst,'linewidth',1.5)

hold on
loglog(xDeg,xlom,'linewidth',1,'color',[0, 0.4470, 0.7410])
%semilogy(xDeg,xPar,'linewidth',1.5)
%plot(xDeg,xPar,'linewidth',1.5)


hold on
loglog(xDeg,xver1,'linewidth',1.2,'color',[0.8500,0.3250, 0.0980])
%semilogy(xDeg,xPow,'linewidth',1.5)
%plot(xDeg,xPow,'linewidth',1.5)


hold on
loglog(xDeg,xver2,'linewidth',1.4,'color',[0.9290, 0.6940, 0.1250])
%semilogy(xDeg,xLog,'linewidth',1.5)
%plot(xDeg,xLog,'linewidth',1.5)

hold on
loglog(xDeg,xver3,'linewidth',1.6,'color',[0, 0.5, 0])
%semilogy(xDeg,xPoC,'linewidth',1.5)
%plot(xDeg,xPoC,'linewidth',1.5)

hold on
loglog(xDeg,xver4,'linewidth',1.8,'color',[1, 0, 0])
%semilogy(xDeg,xExp,'linewidth',1.5)
%plot(xDeg,xExp,'linewidth',1.5)


ylim([0.3 100000]);
xlim([1 10000]);
set(gca,'fontweight','bold','fontsize',12);
xlabel('Node Degree');
ylabel('Frequency');

L=legend('ego-Twitter Network','Pareto Type-I','Power law','Log-normal','Power law cutoff','Exponential','Lomax','GLM Type-I','GLM Type-II','GLM Type-III','GLM Type-IV','Location',[0.5, 0.5, .25, .25]);

%%%%%%%%%%%%%%%%%%%%%%%%%%  Plot of PDFs of all version (GLM Type-I, GLM Type-II,GLM Type-III, GLM Typ-IV) %%%%%%%%%%%%%%%%%%%%%%%%%% 


clear all;

figure

%%%%%%%%%%%%%%%%%%%%%% Plot PDF (GLM TYPE-I) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 1:6;

alpha=1;
beta=0;
gama=1;

term1=(alpha/gama).*((1+x./gama).^(-alpha-1));
term2=(1+(beta./((1+log(1+x./gama)).^2)));
term3=exp(-alpha.*beta.*(log(1+x./gama)./(1+log(1+x./gama))));

y = term1.*term2.*term3;
%plot(x,y,'bp')

xi = linspace(min(x), max(x), 150);                    
yi = interp1(x, y, xi, 'spline', 'extrap');
plot(xi, yi, '-r')

hold on

alpha1=1;
beta1=-0.25;
gama1=1;

term11=(alpha1/gama1).*((1+x./gama1).^(-alpha1-1));
term21=(1+(beta1./((1+log(1+x./gama1)).^2)));
term31=exp(-alpha1.*beta1.*(log(1+x./gama1)./(1+log(1+x./gama1))));

y1 = term11.*term21.*term31;
%plot(x,y1)

xi = linspace(min(x), max(x), 150);                    
yi1 = interp1(x, y1, xi, 'spline', 'extrap');
plot(xi, yi1, '-g')

hold on

alpha2=1;
beta2=-0.5;
gama2=1;

term12=(alpha2/gama2).*((1+x./gama2).^(-alpha2-1));
term22=(1+(beta2./((1+log(1+x./gama2)).^2)));
term32=exp(-alpha2.*beta2.*(log(1+x./gama2)./(1+log(1+x./gama2))));

y2 = term12.*term22.*term32;
%plot(x,y2)

xi = linspace(min(x), max(x), 150);                    
yi2 = interp1(x, y2, xi, 'spline', 'extrap');
plot(xi, yi2, '-b')

hold on

alpha3=1;
beta3=-0.75;
gama3=1;

term13=(alpha3/gama3).*((1+x./gama3).^(-alpha3-1));
term23=(1+(beta3./((1+log(1+x./gama3)).^2)));
term33=exp(-alpha3.*beta3.*(log(1+x./gama3)./(1+log(1+x./gama3))));

y3 = term13.*term23.*term33;
%plot(x,y3)

xi = linspace(min(x), max(x), 150);                    
yi3 = interp1(x, y3, xi, 'spline', 'extrap');
plot(xi, yi3, '-c')

hold on

alpha4=1;
beta4=-0.95;
gama4=1;

term14=(alpha4/gama4).*((1+x./gama4).^(-alpha4-1));
term24=(1+(beta4./((1+log(1+x./gama4)).^2)));
term34=exp(-alpha4.*beta4.*(log(1+x./gama4)./(1+log(1+x./gama4))));

y4 = term14.*term24.*term34;
%plot(x,y4)

xi = linspace(min(x), max(x), 150);                    
yi4 = interp1(x, y4, xi, 'spline', 'extrap');
plot(xi, yi4, '-black')

 ylim([0 0.28]);
 %xlim([0 5]);
 set(gca,'fontsize',14);
 xlabel('x');
 ylabel('f(x)');

 L=legend('\alpha = 1, \gamma = 1, \beta = 0','\alpha = 1, \gamma = 1, \beta = - 0.25','\alpha = 1, \gamma = 1, \beta = - 0.50','\alpha = 1, \gamma = 1, \beta = - 0.75','\alpha = 1, \gamma = 1, \beta = - 0.95','Location',[0.6, 0.6, .25, .25]);

%L.FontWeight='bold';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

figure

%%%%%%%%%%%%%%%%%%%%%% Plot PDF (GLM TYPE-II) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 1:5;

alpha=1;
beta=0;
gama=1;

term1=(1+(beta./(1+(x./gama))));
term2=(alpha./(x+gama));
term3=(1+(x./gama)).^(-alpha);
term4=exp(-alpha.*beta.*((x./gama)./(1+(x./gama))));
  

y = term1.*term2.*term3.*term4;
%plot(x,y,'bp')

xi = linspace(min(x), max(x), 150);                    
yi = interp1(x, y, xi, 'spline', 'extrap');
plot(xi, yi, '-r')

hold on

alpha1=1;
beta1=-0.25;
gama1=1;

term11=(1+(beta1./(1+(x./gama1))));
term21=(alpha1./(x+gama1));
term31=(1+(x./gama1)).^(-alpha1);
term41=exp(-alpha1.*beta1.*((x./gama1)./(1+(x./gama1))));

y1 = term11.*term21.*term31.*term41;
%plot(x,y1)

xi = linspace(min(x), max(x), 150);                    
yi1 = interp1(x, y1, xi, 'spline', 'extrap');
plot(xi, yi1, '-g')

hold on

alpha2=1;
beta2=-0.5;
gama2=1;

term12=(1+(beta2./(1+(x./gama2))));
term22=(alpha2./(x+gama2));
term32=(1+(x./gama2)).^(-alpha2);
term42=exp(-alpha2.*beta2.*((x./gama2)./(1+(x./gama2))));

y2 = term12.*term22.*term32.*term42;
%plot(x,y2)

xi = linspace(min(x), max(x), 150);                    
yi2 = interp1(x, y2, xi, 'spline', 'extrap');
plot(xi, yi2, '-b')

hold on

alpha3=1;
beta3=-0.75;
gama3=1;

term13=(1+(beta3./(1+(x./gama3))));
term23=(alpha3./(x+gama3));
term33=(1+(x./gama3)).^(-alpha3);
term43=exp(-alpha3.*beta3.*((x./gama3)./(1+(x./gama3))));

y3 = term13.*term23.*term33.*term43;
%plot(x,y3)

xi = linspace(min(x), max(x), 150);                    
yi3 = interp1(x, y3, xi, 'spline', 'extrap');
plot(xi, yi3, '-c')

hold on

alpha4=1;
beta4=-0.95;
gama4=1;

term14=(1+(beta4./(1+(x./gama4))));
term24=(alpha4./(x+gama4));
term34=(1+(x./gama4)).^(-alpha4);
term44=exp(-alpha4.*beta4.*((x./gama4)./(1+(x./gama4))));

y4 = term14.*term24.*term34.*term44;
%plot(x,y4)

xi = linspace(min(x), max(x), 150);                    
yi4 = interp1(x, y4, xi, 'spline', 'extrap');
plot(xi, yi4, '-black')

 ylim([0 0.28]);
 %xlim([0 5]);
 set(gca,'fontsize',14);
 xlabel('x');
 ylabel('f(x)');
% 

L=legend('\alpha = 1, \gamma = 1, \beta = 0','\alpha = 1, \gamma = 1, \beta = - 0.25','\alpha = 1, \gamma = 1, \beta = - 0.50','\alpha = 1, \gamma = 1, \beta = - 0.75','\alpha = 1, \gamma = 1, \beta = - 0.95','Location',[0.6, 0.6, .25, .25]);

%L.FontWeight='bold';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

figure

%%%%%%%%%%%%%%%%%%%%%% Plot PDF (GLM TYPE-III) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 1:5;

alpha=1;
beta=0;
gama=1;


term1=(1+(beta.*(log(1+x./gama)))/(x./gama));
term2=((x./gama)./(1+x./gama)).^(beta);
term3=(alpha./(x+gama));
term4=exp(-alpha.*log((1+x./gama).*(((x./gama)/(1+x./gama)).^beta)));
  
y = term1.*term2.*term3.*term4;
%plot(x,y,'bp')

xi = linspace(min(x), max(x), 150);                    
yi = interp1(x, y, xi, 'spline', 'extrap');
plot(xi, yi, '-r')

hold on

alpha1=1;
beta1=-0.25;
gama1=1;

term11=(1+(beta1.*(log(1+x./gama1)))/(x./gama1));
term21=((x./gama1)./(1+x./gama1)).^(beta1);
term31=(alpha1./(x+gama1));
term41=exp(-alpha1.*log((1+x./gama1).*(((x./gama1)/(1+x./gama1)).^beta1)));

y1 = term11.*term21.*term31.*term41;
%plot(x,y1)

xi = linspace(min(x), max(x), 150);                    
yi1 = interp1(x, y1, xi, 'spline', 'extrap');
plot(xi, yi1, '-g')

hold on

alpha2=1;
beta2=-0.5;
gama2=1;

term12=(1+(beta2.*(log(1+x./gama2)))/(x./gama2));
term22=((x./gama2)./(1+x./gama2)).^(beta2);
term32=(alpha2./(x+gama2));
term42=exp(-alpha2.*log((1+x./gama2).*(((x./gama2)/(1+x./gama2)).^beta2)));

y2 = term12.*term22.*term32.*term42;
%plot(x,y2)

xi = linspace(min(x), max(x), 150);                    
yi2 = interp1(x, y2, xi, 'spline', 'extrap');
plot(xi, yi2, '-b')

hold on

alpha3=1;
beta3=-0.75;
gama3=1;


term13=(1+(beta3.*(log(1+x./gama3)))/(x./gama3));
term23=((x./gama3)./(1+x./gama3)).^(beta3);
term33=(alpha3./(x+gama3));
term43=exp(-alpha3.*log((1+x./gama3).*(((x./gama3)/(1+x./gama3)).^beta3)));

y3 = term13.*term23.*term33.*term43;
%plot(x,y3)

xi = linspace(min(x), max(x), 150);                    
yi3 = interp1(x, y3, xi, 'spline', 'extrap');
plot(xi, yi3, '-c')

hold on

alpha4=1;
beta4=-0.95;
gama4=1;


term14=(1+(beta4.*(log(1+x./gama4)))/(x./gama4));
term24=((x./gama4)./(1+x./gama4)).^(beta4);
term34=(alpha4./(x+gama4));
term44=exp(-alpha4.*log((1+x./gama4).*(((x./gama4)/(1+x./gama4)).^beta4)));

y4 = term14.*term24.*term34.*term44;
%plot(x,y4)

xi = linspace(min(x), max(x), 150);                    
yi4 = interp1(x, y4, xi, 'spline', 'extrap');
plot(xi, yi4, '-black')

 ylim([0 0.28]);
 %xlim([0 5]);
 set(gca,'fontsize',14);
 xlabel('x');
 ylabel('f(x)');
% 
L=legend('\alpha = 1, \gamma = 1, \beta = 0','\alpha = 1, \gamma = 1, \beta = - 0.25','\alpha = 1, \gamma = 1, \beta = - 0.50','\alpha = 1, \gamma = 1, \beta = - 0.75','\alpha = 1, \gamma = 1, \beta = - 0.95','Location',[0.6, 0.6, .25, .25]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

figure

%%%%%%%%%%%%%%%%%%%%%% Plot PDF (GLM TYPE-IV) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = 1:5;

alpha=1;
beta=0;
sigma=1;

term1 = (alpha/sigma);
term2 = (beta+1+log(1+x./sigma))./(1+(x./sigma));
term3 = ((log(1+x./sigma)).^beta)./((1+(log(1+x./sigma))).^(beta+1));
term4 = exp(-alpha*(((log((x./sigma)+1)).^(beta+1))./((log((x./sigma)+1)+1).^(beta))));
y = term1.*term2.*term3.*term4;
%plot(x,y,'bp')

xi = linspace(min(x), max(x), 150);                    
yi = interp1(x, y, xi, 'spline', 'extrap');
plot(xi, yi, '-r')

hold on

alpha1=1;
beta1=-0.25;
sigma1=1;

term11 = (alpha1/sigma1);
term21 = (beta1+1+log(1+x./sigma1))./(1+(x./sigma1));
term31 = ((log(1+x./sigma1)).^beta1)./((1+(log(1+x./sigma1))).^(beta1+1));
term41 = exp(-alpha1*(((log((x./sigma1)+1)).^(beta1+1))./((log((x./sigma1)+1)+1).^(beta1))));
y1 = term11.*term21.*term31.*term41;
%plot(x,y1)

xi = linspace(min(x), max(x), 150);                    
yi1 = interp1(x, y1, xi, 'spline', 'extrap');
plot(xi, yi1, '-g')

hold on

alpha2=1;
beta2=-0.5;
sigma2=1;

term12 = (alpha2/sigma2);
term22 = (beta2+1+log(1+x./sigma2))./(1+(x./sigma2));
term32 = ((log(1+x./sigma2)).^beta2)./((1+(log(1+x./sigma2))).^(beta2+1));
term42 = exp(-alpha2*(((log((x./sigma2)+1)).^(beta2+1))./((log((x./sigma2)+1)+1).^(beta2))));
y2 = term12.*term22.*term32.*term42;
%plot(x,y2)

xi = linspace(min(x), max(x), 150);                    
yi2 = interp1(x, y2, xi, 'spline', 'extrap');
plot(xi, yi2, '-b')

hold on

alpha3=1;
beta3=-0.75;
sigma3=1;

term13 = (alpha3/sigma3);
term23 = (beta3+1+log(1+x./sigma3))./(1+(x./sigma3));
term33 = ((log(1+x./sigma3)).^beta3)./((1+(log(1+x./sigma3))).^(beta3+1));
term43 = exp(-alpha3*(((log((x./sigma3)+1)).^(beta3+1))./((log((x./sigma3)+1)+1).^(beta3))));
y3 = term13.*term23.*term33.*term43;
%plot(x,y3)

xi = linspace(min(x), max(x), 150);                    
yi3 = interp1(x, y3, xi, 'spline', 'extrap');
plot(xi, yi3, '-c')

hold on

alpha4=1;
beta4=-0.95;
sigma4=1;

term14 = (alpha4/sigma4);
term24 = (beta4+1+log(1+x./sigma4))./(1+(x./sigma4));
term34 = ((log(1+x./sigma4)).^beta4)./((1+(log(1+x./sigma4))).^(beta4+1));
term44 = exp(-alpha4*(((log((x./sigma4)+1)).^(beta4+1))./((log((x./sigma4)+1)+1).^(beta4))));
y4 = term14.*term24.*term34.*term44;
%plot(x,y4)

xi = linspace(min(x), max(x), 150);                    
yi4 = interp1(x, y4, xi, 'spline', 'extrap');
plot(xi, yi4, '-black')

 ylim([0 0.28]);
 %xlim([0 5]);
 set(gca,'fontsize',14);
 xlabel('x');
 ylabel('f(x)');
% 
L=legend('\alpha = 1, \gamma = 1, \beta = 0','\alpha = 1, \gamma = 1, \beta = - 0.25','\alpha = 1, \gamma = 1, \beta = - 0.50','\alpha = 1, \gamma = 1, \beta = - 0.75','\alpha = 1, \gamma = 1, \beta = - 0.95','Location',[0.6, 0.6, .25, .25]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
