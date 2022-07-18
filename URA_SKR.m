clc;close all;clear all

Mx= 6;      % Element num on x-axis
My= 6;      % Element num on x-axis
M= Mx*My;   % Element num

% Parameters
fc = 30e9;
c = 3e8;
lam = c/fc; % Wavelength
d=lam/2;    % Element spacing
POP_SIZE = 100000;

% Alice and Bob angles
in1phi=40;in1the=70;
out1phi=120;out1the=60;

in2phi=120;in2the=60;
out2phi=40;out2the=70;

% Array settings
n=Mx;
xPos=linspace(0,(n-1)*d,n);
yPos=linspace(0,(n-1)*d,n);
[YPOS,XPOS]=meshgrid(xPos,yPos);
XPOS=reshape(XPOS,[1,numel(XPOS)]);
YPOS=reshape(YPOS,[1,numel(YPOS)]);
ZPOS = zeros(1, M);


%% Pattern and mutual information
SNR=0:5:20;    % SNR range
MI = zeros(1,length(SNR));MI_th = zeros(1,length(SNR));
MI_RSS=zeros(1,length(SNR));MI_phase=zeros(1,length(SNR));
for m=1:length(SNR)
    % Weights
    pop=rand(POP_SIZE,M)*2*pi;
    w = exp(1j*pop);
    
    % Patterns
    Pattern1=zeros(POP_SIZE,1); % Alice to Bob
    Pattern2=zeros(POP_SIZE,1); % Bob to Alice
    for i=1:POP_SIZE
        weight= w(i,:);
        % Alice to Bob
        Pattern1(i,:) = NewarrayFactor(XPOS, YPOS, ZPOS, weight, fc, c, out1the, out1phi,0,0,in1the,in1phi);
        % Bob to Alice
        Pattern2(i,:) = NewarrayFactor(XPOS, YPOS, ZPOS, weight, fc, c, out2the, out2phi,0,0,in2the,in2phi);
    end
    
    % Add noise
    noise_power = 1/power(10,SNR(m)/10)/2;
    noise_real= randn(POP_SIZE,1)*sqrt(noise_power)+1i*randn(POP_SIZE,1)*sqrt(noise_power);
    Pattern1 = Pattern1 + noise_real;
    noise_imag= randn(POP_SIZE,1)*sqrt(noise_power)+1i*randn(POP_SIZE,1)*sqrt(noise_power);
    Pattern2 = Pattern2 + noise_imag;
    
    % Mutual Information using ITE toolbox
    ds = [1;1];
    mult = 1;
    % Real part MI
    Y_real = [transpose(real(Pattern1));transpose(real(Pattern2))];
    co = IShannon_AP2_initialization(mult);
    MI_real = IShannon_AP_estimation(Y_real,ds,co)*log2(exp(1));
    MI_RSS(1,m)=MI_real;
    % Imag part MI
    Y_imag = [transpose(imag(Pattern1));transpose(imag(Pattern2))];
    co = IShannon_AP_initialization(mult);
    MI_imag = IShannon_AP2_estimation(Y_imag,ds,co)*log2(exp(1));
    MI_phase(1,m)=MI_imag;
    
    MI(1,m) = MI_real+MI_imag;
    
    % Theoretical MI
    MI_th(1,m) = log2(1+M/2/(noise_power*2+noise_power*noise_power/(M/2)));
end

%% Plot
plot(SNR,MI,'-o')
hold on
plot(SNR,MI_th,'-x')
xlabel('Signal-to-noise Ratio (dB)')
ylabel('Secret Key Rate')
legend('Simulation Result','Analytical Result')
grid on
ylim([0 16])
grid on