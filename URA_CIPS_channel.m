clc;close all;clear all

%% URA RIS settings
Mx= 6;      % Element number on x-axis
My= 6;      % Element number on y-axis
M= Mx*My;   % Total number of elements

fc = 30e9;  % Working frequency
c = 3e8;
lam = c/fc; % Wavelength
d=lam/2;    % Element spacing
POP_SIZE = 100000; % Random/Monte Carlo weight size

% Alice and Bob locating angles with respect to RIS (can be changed)
in1phi=50;      in1the=60;
out1phi=150;    out1the=80;
in2phi=100;     in2the=40;
out2phi=280;    out2the=10;

%% Generate patterns
% random weights
pop=rand(POP_SIZE,M)*2*pi;
w = exp(1j*pop);
n=Mx;

% Element positions
xPos=linspace((-n/2+0.5)*d,(n/2-0.5)*d,n); 
yPos=linspace((-n/2+0.5)*d,(n/2-0.5)*d,n);
[YPOS,XPOS]=meshgrid(xPos,yPos);
XPOS=reshape(XPOS,[1,numel(XPOS)]);  
YPOS=reshape(YPOS,[1,numel(YPOS)]);
ZPOS = zeros(1, M); 

% Patterns
Pattern1=zeros(POP_SIZE,1);
Pattern2=zeros(POP_SIZE,1);
for i=1:POP_SIZE
    weight= w(i,:);
    Pattern1(i,:) = NewarrayFactor(XPOS, YPOS, ZPOS, weight, fc, c, out1the, out1phi,0,0,in1the,in1phi);
    Pattern2(i,:) = NewarrayFactor(XPOS, YPOS, ZPOS, weight, fc, c, out2the, out2phi,0,0,in2the,in2phi);
end

%% Plot magnitude p.d.f
pts_abs = (min(abs(Pattern1)):0.5:max(abs(Pattern1)));
[f_abs1,xi_abs1] = ksdensity(abs(Pattern1),pts_abs);
[f_abs2,xi_abs2] = ksdensity(abs(Pattern2),pts_abs);

% Theoretical p.d.f
sig = M/2;
y = xi_abs1/sig.*exp(-power(xi_abs1,2) / (2 *sig)) ;
figure
line1=plot(xi_abs1,f_abs1,'-- ^');
line1.MarkerIndices = 1:4:length(xi_abs1);

hold on 
line2=plot(xi_abs2,f_abs2,'-- d');
line2.MarkerIndices = 1:3:length(xi_abs1);

hold on 
line3=plot(xi_abs1,y,'-- x');
line3.MarkerIndices = 1:2:length(xi_abs1);

xlabel('Magnitude')
ylabel('Probability Density')
legend('Pattern 1','Pattern 2','Analytical result')
grid on

%% Plot phase p.d.f
angle1=angle(Pattern1);
angle2=angle(Pattern2);

a=find(angle1<0);
angle1(a)=angle1(a)+pi*2;
b=find(angle2<0);
angle2(b)=angle2(b)+pi*2;

pts_phase = 0:0.05:2*pi;
[f_phase1,xi_phase1] = ksdensity(angle1,pts_phase,'Bandwidth',0.03);
[f_phase2,xi_phase2] = ksdensity(angle2,pts_phase,'Bandwidth',0.03);
figure
line4=plot(xi_phase1,f_phase1,'-- ^');
line4.MarkerIndices = 5:10:length(xi_phase1);

hold on
line5=plot(xi_phase2,f_phase2,'-- d');
line5.MarkerIndices = 5:12:length(xi_phase1);

hold on
line=1/2/pi*ones(1,length(f_phase1));
line6=plot(xi_phase1,line,'-- x');
line6.MarkerIndices = 5:8:length(xi_phase1);

set(gca,'Xtick',[0,pi/2,pi,3*pi/2,2*pi]);
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
xlim([0 2*pi])
ylim([0 1/2/pi+0.02])
xlabel('Phase (arc)')
ylabel('Probability Density')
legend('Pattern 1','Pattern 2','Analytical result')
grid on

