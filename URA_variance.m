clc;close all;clear all

%% Parameters
Mx= 6:2:16;%Element quantity
My= Mx;
M= Mx.*My;

fc = 30e9;
c = 3e8;
lam = c/fc; % Wavelength
d=lam/2;	% Element spacing
POP_SIZE = 10000;

%  Alice and Bob angles with respect to RIS
in1phi=40;in1the=30;
out1phi=70;out1the=100;
in2phi=110;in2the=20;
out2phi=310;out2the=10;

%% Variance
variance1=zeros(length(M),1);
variance2=zeros(length(M),1);
for M_num=1:length(M)
    m=M(M_num);
    
    % Weights
    pop=rand(POP_SIZE,m)*2*pi;
    w = exp(1j*pop);
    
    % Array element positions
    n=Mx(M_num);
    xPos=linspace((-n/2+0.5)*d,(n/2-0.5)*d,n); 
    yPos=linspace((-n/2+0.5)*d,(n/2-0.5)*d,n);
    [YPOS,XPOS]=meshgrid(xPos,yPos);
    XPOS=reshape(XPOS,[1,numel(XPOS)]);  
    YPOS=reshape(YPOS,[1,numel(YPOS)]);
    ZPOS = zeros(1, m); 
    
    % Patterns
    Pattern1=zeros(POP_SIZE,1);
    Pattern2=zeros(POP_SIZE,1);
    for i=1:POP_SIZE
        weight= w(i,:);
        Pattern1(i,:) = NewarrayFactor(XPOS, YPOS, ZPOS, weight, fc, c, out1the, out1phi,0,0,in1the,in1phi);
        Pattern2(i,:) = NewarrayFactor(XPOS, YPOS, ZPOS, weight, fc, c, out2the, out2phi,0,0,in2the,in2phi);
    end
    variance1(M_num)=var(Pattern1);
    variance2(M_num)=var(Pattern2);
end

%% Plot
figure
plot(M,variance1,'-- ^')
xlabel('Number of Elements')
ylabel('Channel Variance')
grid on
hold on
plot(M,variance2,'-- >')

% Theoretical analysis
hold on
line = plot(M,M,'-- o');
legend('Angle pair 1','Angle pair 2','Analytical result')