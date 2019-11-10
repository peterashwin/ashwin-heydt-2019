
%%
% PA & AvdH
% Matlab script for
% ECS + tipping
%
% Oct 2019
%
%
%% Set seed
clear variables;
rng(1811);

% output fig names
fname=sprintf('ECST-nov-2019-%s',date);

printfigs=true;
savedata=false;

% simulation time
total=50000;

%% params for functions
% temp thresholds (K)
p.T1=278;
p.T2=288;
% albedo for low and high temp limit
p.AB1=0.52;
p.AB2=0.47;
% noise for TG and for ln(C)
p.etaTG=5e-6;
p.etaLC=2e-6;
% 1/EPSH is steepness of switch
p.EPSH=3.0; %0.273 K
% QG is solar input
p.QG=342.0;

%% CO2 atmospheric coeffs
%p.G=1.35e2;
p.C0=280;
% CO2 atmospheric processes for low and high temp limits
p.AT=5.35;
% emissivity switch
p.EM1=0.53;
p.EM2=0.39;
p.EPSE=20.0;
p.T0=288;

%% cutoff for random wandering of CO2 conc
p.CMIN=100;
p.CMAX=1000;
p.CK=1e-7;

%% thermal interia
p.CT=5e8;

%% Boltzmann coeff
p.sigmaB=5.67e-8;

%% Scaling for time units (1/sec)
% years per sec
p.tau=1/(60*60*24*365);

%% timestep
st=5;
h=1/st;
ni=st*total;

%% dynamic variables with initial conditions
% initial temp (K), log atmospheric CO2 (ppm)

yy=[290 log(300)];

t=0;
i=1;
y=zeros(ni,2);
tt=zeros(ni,1);
y(i,:)=yy;
tt(i)=t;

while (t<total)
    % Euler-Maruyama step
    fftemp=FF(yy,p);
    ynew=yy+h*fftemp/p.tau+...
        sqrt(h/p.tau)*[p.etaTG p.etaLC].*randn([1 2]);
    t=t+h;
    i=i+1;
    y(i,:)=ynew;
    tt(i)=t;
    yy=ynew;
end
itot=i-1;
% global mean temp
TG=y(:,1);
% atmos CO2
CO=exp(y(:,2));

fignum=1;
%fignum=fignum+1;
%% figs of alpha, epsilon and Gamma
f1=figure(fignum);
clf;
f1.PaperUnits='centimeters';
f1.PaperSize=[12 12];
f1.Units='centimeters';
f1.InnerPosition=[1 1 12 12];

clf;
Tmin=270; %Tmin=240;
Tmax=300; %Tmax=330;

%transient=100;
%plot(CO(transient:end),TG(transient:end)-273.15);
% plot equilibrium, also alpha and epsilon
ts=[Tmin:0.1:Tmax];
tsc=ts-273.15;
DR=zeros(size(ts));
CC=zeros(size(ts));
al=zeros(size(ts));
ep=zeros(size(ts));
for ii=1:size(ts,2)
    DR(ii)=RCOCurve(ts(ii),p);
    CC(ii)=p.C0*exp(DR(ii)/p.AT);
    al(ii)=ALPHA(ts(ii),p);
    ep(ii)=EPSILON(ts(ii),p);
end

subplot(2,2,1)
plot(tsc,al);
hold on;
xlabel('T [C]');
ylabel('\alpha(T)')
xlim([tsc(1) tsc(end)])
ylim([0.45 0.55]);
text(0.9,0.05,'a','Units','normalized');
hold off;

subplot(2,2,2)
plot(tsc,ep);
hold on;
xlabel('T [C]');
ylabel('\epsilon(T)')
xlim([tsc(1) tsc(end)])
ylim([0.35 0.55]);
text(0.9,0.05,'b','Units','normalized');
hold off;

subplot(2,2,3)
plot(CC,tsc);
hold on;
xlabel('C [ppm]')
ylabel('T [C]');
ylim([tsc(1) tsc(end)])
xlim([0 1e3]);
text(0.9,0.05,'c','Units','normalized');
hold off;

subplot(2,2,4)
plot(DR,tsc);
hold on;
xlabel('\Delta R_{CO_2} [W/(m^2/s)]')
ylabel('T [C]');
ylim([tsc(1) tsc(end)])
xlim([-10 10]);
%xlim([0 1.1*p.CMAX]);
text(0.9,0.05,'d','Units','normalized');
hold off;

savefigure(fname,fignum,printfigs);

%keyboard


%% figure of timeseries
figure(fignum);
clf;
%
subplot(3,1,1)
plot(tt,TG-273.15);
xlim([0, total]);
xlabel('t [yr]')
ylabel('T(t) [C]');
text(0.95,0.1,'a','Units','normalized');
%
subplot(3,1,2)
plot(tt,CO);
hold on;
plot(tt,p.CMAX*ones(size(tt)),'r');
plot(tt,p.CMIN*ones(size(tt)),'r');
xlim([0, total]);
ylim([0 1600]);
xlabel('t [yr]')
ylabel('C(t) [ppm]');
text(0.95,0.1,'b','Units','normalized');
hold off;
%
subplot(3,1,3)
plot(tt,log(CO/p.C0));
hold on;
plot(tt,log(p.CMAX/p.C0)*ones(size(tt)),'r');
plot(tt,log(p.CMIN/p.C0)*ones(size(tt)),'r');
xlim([0, total]);
ylim([-2 2]);
xlabel('t [yr]')
ylabel('ln(C(t)/C_0)');
text(0.95,0.1,'c','Units','normalized');
hold off;

savefigure(fname,fignum,printfigs);

fignum=fignum+1;
%% figure of T vs C
f2=figure(fignum);
clf;
f2.PaperUnits='centimeters';
f2.PaperSize=[12 7];
f2.Units='centimeters';
f2.InnerPosition=[1 1 12 7];
transient=100;
subplot(1,2,1);
plot(CO(transient:end),TG(transient:end)-273.15);
xlabel('C [ppm]')
ylabel('T(t) [C]');
text(0.9,0.1,'a','Units','normalized');
hold on;
% plot null cline
ts=[min(y(transient:end,1))-2:0.1:max(y(transient:end,1))+2];
tsc=ts-273.15;
CC=zeros(size(ts));
for ii=1:size(ts,2)
    CC(ii)=p.C0*exp((RCOCurve(ts(ii),p))/p.AT);
end
plot(CC,tsc);
ylim([tsc(1) tsc(end)])
xlim([0 1.1*p.CMAX]);
hold off;

%savefigure(fname,fignum,printfigs);

%fignum=fignum+1;
%% figure of T vs Delta R 
%figure(fignum);
%clf;
subplot(1,2,2);
transient=100;
RCO2=p.AT*log(exp(y(:,2))/p.C0);
plot(RCO2(transient:end),TG(transient:end)-273.15);
xlabel('\Delta R_{[CO_2]}(t) [W/m^2]')
ylabel('T(t) [C]');
hold on;
% plot null cline
ts=min(y(transient:end,1))-2:0.2:max(y(transient:end,1))+2;
tsc=ts-273.15;
RCurve=zeros(size(ts));
for ii=1:size(ts,2)
    RCurve(ii)=RCOCurve(ts(ii),p);
end
plot(RCurve,tsc);
ylim([tsc(1) tsc(end)])
xlim([-10 10]);
text(0.9,0.1,'b','Units','normalized');
hold off;

savefigure(fname,fignum,printfigs);

%% downsample to one per time unit
ttt=tt(1:st:end);
TGT=TG(1:st:end);
RCOT=RCO2(1:st:end);

fignum=fignum+1;
%% figs of local sensitivities etc
f1=figure(fignum);
clf;
f1.PaperUnits='centimeters';
f1.PaperSize=[10 15];
f1.Units='centimeters';
f1.InnerPosition=[1 1 10 15];

TCC=zeros(size(ttt));
SS=zeros(size(ttt));
for ii=1:length(ttt)
    % find temp on equilibrium curve with this RCO2
     temp=TCurve(TGT(ii),RCOT(ii),p);
     TCC(ii)=temp;
     % find dRCO2/dT 
     temp2=dRCOCurve(temp,p);
     % inverse is sensitivity
     SS(ii)=1/temp2;
end
%%
p1=subplot(5,1,1);
plot(ttt,RCOT);
hold on;
plot([0 total],-1.744*[1 1],'r');
plot([0 total],3.004*[1 1],'r');
hold off;
xlim([0, total]);
%xlabel('t [yr]')
ylabel('\Delta R');
text(0.95,0.8,'a','Units','normalized');
p1.XTickLabel=[];
%
%%
p2=subplot(5,1,2);
plot(ttt,TGT-273.15);
xlim([0, total]);
ylim([-5 25]);
%xlabel('t [yr]')
ylabel('T');
p2.XTickLabel=[];
text(0.95,0.8,'b','Units','normalized');
%
p3=subplot(5,1,3);
plot(ttt,SS);
xlim([0, total]);
ylim([0 2.0]);
%xlabel('t [yr]')
ylabel('S');
%p3.YTick=([0 0.5 1.0 1.5 2.0]);
p3.XTickLabel=[];
text(0.95,0.8,'c','Units','normalized');

%
% size of window
nw=500;


np=length(ttt);
ttc=ttt(nw:np);
% moving average T
TGA=zeros(size(ttt));
for i=nw:np
    TGClip=TGT(i+1-nw:i);
    RCOClip=RCO2(i+1-nw:i);
    % moving average
    TGA(i)=mean(TGClip);
end

% moving SD
TGSD=zeros(size(ttt));
% moving AR
TGAR1=zeros(size(ttt));

% delay for AR
nd=1;
for n=nw:np-nd
    % moving average
    temp=sum((TGT(n-nw+1:n)-TGA(n)).^2);
    % SD
    TGSD(n)=sqrt(temp/nw);
    % AR(nd)
    TGAR1(n)=sum((TGT(n-nw+nd:n-1)-TGA(n)).*(TGT(n-nw+1+nd:n)-TGA(n)))/temp;
end     

%%
p4=subplot(5,1,4);
plot(ttt,TGSD);
hold on;
xlim([0, total]);
ylim([-2 15]);
%xlabel('t [yr]')
ylabel('sd(t)');
p4.XTickLabel=[];
text(0.95,0.8,'d','Units','normalized');

hold off;
%
%%
subplot(5,1,5)
plot(ttt,TGAR1);
hold on;
xlim([0, total]);
ylim([0.97 1.005]);
xlabel('t [yr]')
ylabel('AR1(t)');
text(0.95,0.8,'e','Units','normalized');
%
hold off;

savefigure(fname,fignum,printfigs);


%%
if(savedata==true)
    % sample data one per tme unit
    sel=1:st:itot;
    rundata=[tt(sel),TG(sel),CO(sel),RCO2(sel)];
    save(sprintf('./data/data-%s.txt',fname),'rundata','-ascii');
end

keyboard


function temp=H(x,p)
%% smoothed heaviside
temp=(1+tanh(x./p.EPSH))./2.0;
end

function temp=ALPHA(x,p)
%% switch albedo
temp=p.AB1.*H(p.T1-x,p)+p.AB2.*H(x-p.T2,p)+...
    (p.AB1+(p.AB2-p.AB1)*(x-p.T1)./(p.T2-p.T1)).*H(x-p.T1,p).*H(p.T2-x,p);
end

function temp=EPSILON(x,p)
%% switch EM
tt=(1+tanh((x-p.T0)/p.EPSE))/2;
temp=p.EM1+(p.EM2-p.EM1).*tt;
end

function dy=FF(y,p)
%% RHS of ODE
% Energy balance from Dijkstra and Viebahn 2015
% with tanh emissivity variation
%
TG=y(1);
LC=y(2);

tt=p.QG*(1-ALPHA(TG,p))+p.AT*log(exp(LC)/p.C0)-...
    p.sigmaB*EPSILON(TG,p)*TG^4;
dy(1) = (tt)/p.CT;
% random walk of CO2 with soft constraints
dy(2) = p.CK*F(LC,p);
end

function temp=F(x,p)
% soft constraints on random walk
LC1=log(p.CMIN);
LC2=log(p.CMAX);
temp=(heaviside(LC1-x).*(LC1-x)+heaviside(x-LC2).*(LC2-x));
end

function temp=RCOCurve(T,p)
% equilibrium RCO2 as function of temp
temp=-p.QG.*(1-ALPHA(T,p))...
        +p.sigmaB.*EPSILON(T,p).*T.^4;
end

function temp=dRCOCurve(T,p)
% d/dT of RCOCurve by finite difference
h=1e-6;
temp=(RCOCurve(T+h,p)-RCOCurve(T,p))/h;
end

function temp=TCurve(Tin,RCO,p)
% equilibrium T as function of RCO2 near Tin
% simple ode search for nearest equilibrium
ts=[0 10];
[tt,yy]=ode15s(@(t,y) RCO-RCOCurve(y,p),ts,Tin);
temp=yy(end);

end
