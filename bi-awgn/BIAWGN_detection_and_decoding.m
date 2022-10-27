function BIAWGN_detection_and_decoding(snr, n, Efa, Emd, Eie)
%
% BIAWGN_detection_and_decoding(rho, n, Efa, Emd, Eie):
% Computation of different detection and decoding bounds for
% the binary-input AWGN channel.
%
% INPUTS:
% rho   = SNR [dB]
% n     = Blocklength
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% If debug=1, this function plots all the bounds with different QY.
% If debug~=1, the function only saves the results of the bounds.

debug = 1;

if debug
    snr = 3;
    n = 10:10:300;
    Efa = 1e-4;
    Emd = 1e-4;
    Eie = 1e-3;
end
% snr = 10*log10(rho);
rho = 10^(snr/10);

% INITIALIZATIONS:
R_ach_th1 = nan(size(n));
R_ach_th1_pow = nan(size(n));
R_conv_Th2 = nan(size(n));
R_ach_pre = zeros(size(n));
R_ach_pre_pow = zeros(size(n));
R_conv_pre = zeros(size(n));
R_conv_pre_pow = zeros(size(n));
R_genie_conv = nan(size(n));
R_genie_ach_DT = nan(size(n));
npre_ach = zeros(size(n));
npre_conv = zeros(size(n));
Px_th1 = nan(size(n));
Px_th2 = nan(size(n));
pow_th1 = nan(size(n));
P1_ach_pre = nan(size(n));
P1_conv_pre = nan(size(n));

disp(['Simulation for n = ' num2str(n)] )

for i = 1:length(n)
    
%     R_imperfect_ach_X(i) = Achievability_X(n(i), rho, Efa, Emd, Eie);
%     disp(['Detection and decoding achievability (Th1 Qy=PY): (R, n) = (' num2str(R_imperfect_ach_X(i)) ',' num2str(n(i)) ')']);
    
    [R_ach_th1(i),Px_th1(i)] = Achievability_Th1(n(i), rho, Efa, Emd, Eie);
    disp(['Detection and decoding achievability (Th1): (R, n, Px) = (' num2str(R_ach_th1(i)) ',' num2str(n(i)) ',' num2str(Px_th1(i)) ')']);
    
%     [R_imperfect_ach_th1_pow(i),pow_th1(i)] = Achievability_Th1_pow(n(i), rho, Efa, Emd, Eie);
%     disp(['Detection and decoding achievability (Th1 with power allocation): (R, n, P1) = (' num2str(R_imperfect_ach_th1_pow(i)) ',' num2str(n(i)) ',' num2str(pow_th1(i)) ')']);
    
    [R_conv_Th2(i),Px_th2(i)] =  converse_Th2(n(i), rho, Efa, Emd, Eie);
    disp(['Detection and decoding converse (Th2 QY=PY): (R, n) = (' num2str(R_conv_Th2(i)) ',' num2str(n(i)) ',' num2str(Px_th2(i)) ')']);
    
    [R_ach_pre(i), npre_ach(i)] = Achievability_Preamble(n(i), rho, Efa, Emd, Eie);
    disp(['Detection and decoding achievability preamble (Th2): (R, n, npre) = (' num2str(R_ach_pre(i)) ',' num2str(n(i)) ',' num2str(npre_ach(i)) ')']);
    
%     [R_imperfect_ach_pre_pow(i),P1_ach_pre(i)] = Achievability_Preamble_pow(n(i), rho, Efa, Emd, Eie);
%     disp(['Detection and decoding achievability preamble (Th2 with power allocation): (R, n, P1) = (' num2str(R_imperfect_ach_pre_pow(i)) ',' num2str(n(i)) ',' num2str(P1_ach_pre(i)) ')']);
    
    [R_conv_pre(i), npre_conv(i)] = Converse_Preamble(n(i), rho, Efa, Emd, Eie);
    disp(['Detection and decoding converse preamble (Th3): (R, n) = (' num2str(R_conv_pre(i)) ',' num2str(n(i)) ',' num2str(npre_conv(i)) ')']);
    
%     [R_imperfect_conv_pre_pow(i),P1_conv_pre(i)] = Converse_Preamble_pow(n(i), rho, Efa, Emd, Eie);
%     disp(['Detection and decoding converse preamble (Th3 with power allocation): (R, n, P1) = (' num2str(R_imperfect_conv_pre_pow(i)) ',' num2str(n(i)) ',' num2str(P1_ach_pre(i)) ')']);
%     
    R_genie_conv(i) = metaconverse_biawgn_fixed_s(Eie,n(i),rho,1);
    disp(['Genie-detection and decoding converse (metaconverse): (R, n) = (' num2str(R_genie_conv(i)) ',' num2str(n(i)) ')']);
%     
    R_genie_ach_DT(i) = perfect_achievability_DT(n(i), rho, Eie);
    disp(['Genie-detection and decoding achievability (DT): (R, n) = (' num2str(R_genie_ach_DT(i)) ',' num2str(n(i)) ')']);
    
end

C = capacity(rho)*ones(size(n));
[n, R_NA] = NA(rho,n,Eie);

if debug == 1
    figure
    plot(n, R_conv_Th2,'b'); hold on;
    plot(n, R_ach_th1,'--b');
%     plot(n, R_ach_th1_pow,'color',[1 0.5 0.5]);
    plot(n, R_ach_pre,'--','color',[0 0.8 0]);hold on;
%     plot(n, R_ach_pre_pow,'color',[0.7 0.8 0.1]);hold on;
    plot(n, R_conv_pre,'color',[0 0.8 0]);hold on;
%     plot(n, R_conv_pre_pow,'color',[0.5 0 0.5]);hold on;
    plot(n, R_genie_conv,'r');
    plot(n, R_genie_ach_DT,'--r');
    % plot(n, C,'k');
    % plot(n,R_NA);
    legend('Det-dec converse (Th2)','Det-dec ach. (Th1)','Det-dec ach. preamble (based Th1)','Det-dec converse preamble (based Th2)','genie metaconverse','genie DT')
    title(['Bi-AWGN Parameters: ed=' num2str(Eie) ', emd=' num2str(Emd) ', efa=' num2str(Efa)])
    ylabel('rate, R [bits]')
    xlabel('blocklnegth, n')
    
    figure
    plot(n,Px_th2); hold on
    plot(n,Px_th1)
    title(['Bi-AWGN Parameters: eie=' num2str(Eie) ', emd=' num2str(Emd) ', efa=' num2str(Efa)])
    ylabel('Px')
    xlabel('blocklnegth, n')
    legend('Imp. det. converse (Th3 QY=PY)','Imp. det. ach. (Th2)')

end

filename = ['BIAWGN_DET_SNR_' num2str(snr) 'dB_ed_' sprintf('%.0e', Eie) '_efa_' sprintf('%.0e', Efa) '_emd_' sprintf('%.0e', Emd) ];
data = [n; R_ach_th1 ; R_ach_th1; R_ach_th1_pow; R_ach_pre; R_ach_pre_pow; R_conv_Th2; R_conv_pre; R_conv_pre_pow; R_genie_conv;R_genie_ach_DT;R_NA;C;npre_ach; npre_conv; Px_th1; Px_th2; pow_th1; P1_ach_pre; P1_conv_pre]';

csvwrite([filename '.csv'],data)
save([filename '.mat'],'data')
end

function [R, Px_1] = Achievability_Th1(n, rho, Efa, Emd, Eie)
%
% Achievability_Th1(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 2 for
% detection and decoding for the bi-AWGN channel. This achievability skew
% the input distribution to facilitate detection. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR [dB]
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% Px_1 = Prob[X=1]

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

% Definitions according to the paper:
i_XY_fcn = @(Px, x, y) log(exp(y.*x) ./ (Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) ;
r_XY_fcn = @(Px, y) -log((Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) + rho/2 ;
% r_XY_fcn2 = @(Px, y) -log((Px * exp(y*sqrt(rho)+rho) + (1-Px) * exp(-y*sqrt(rho)-rho))) + rho/2 ;

% STEP 1: Skew input: try to find Px such that emd and efa are satisfied.
steplength = 0.01;
Px_1 = 0.5-steplength;
emd=1;%initital values
while emd > Emd
    r_Y_noise = 0;
    %     r_Y = 0;
    Px_1 = Px_1 + steplength;
    %generate random variable r(Y). (Definition in the paper)
    for i = 1:n
        Y = gaussian_rv(0,1,1,N);
        r_Y_noise = r_Y_noise + r_XY_fcn(Px_1, Y);
        %r_Y = r_Y + r_XY_fcn2(Px_1, Y);
    end
    
    %Obtain ccdf of the r.v.'s
    [cdf_R, Rx] = ecdf(-r_Y_noise); %use built in function to obtain empirical cdf
    Rx(1) = Rx(1) + 1e-8;
    gamma_fa = interp1(cdf_R, Rx, 1-Efa); % this is the minimum value of the threshold we can use
    
    %[cdf_R2, Rx2] = ecdf(-r_Y); %use built in function to obtain empirical cdf
    %Rx2(1) = Rx2(1) + 1e-8;
    %emd2 = interp1(Rx2,cdf_R2, gamma_fa); % This version requires many more samples
    emd = mean(exp(-r_Y_noise) .* (-r_Y_noise < gamma_fa)); % Includes a change of measure.
    
    if Px_1 >= 1
        R=0; return;
    end
end

% STEP 2: Find the maximum rate achievabable. Extra skewness of the input
% allowed.
R = 0;
R_test = 0;
Px_1 = Px_1 - steplength;
while R <= R_test
    r_Y_noise = 0;
    i_XY = 0;
    Px_1 = Px_1 + steplength;
    if Px_1 > 1
        Px_1 = 1;
        return;
    end
    %generate random variable r(Y): (Definition in the paper)
    for i = 1:n
        %X = sqrt(rho)*(-1).^(rand(1,N) > Px); %Pr[X=1]=Px
        Y = gaussian_rv(0,1,1,N);
        r_Y_noise = r_Y_noise + r_XY_fcn(Px_1, Y);
    end
    
    [cdf_R, Rx] = ecdf(r_Y_noise); %use built in function to obtain empirical cdf
    %Rx(1) = Rx(1) + 1e-8;
    gamma_fa = interp1(cdf_R, Rx, Efa); % this is the minimum value of the threshold we can use
    
    emd = mean(exp(-r_Y_noise) .* (r_Y_noise > gamma_fa));
    
    % Generate information density i(X;Y): (Definition in the paper)
    for i = 1:n
        X = sqrt(rho)*(-1).^(rand(N,1) > Px_1); %Pr[X=1]=Px
        Y = X + gaussian_rv(0,1,1,N);
        i_XY = i_XY + i_XY_fcn(Px_1, X, Y);
    end
    % STEP 3: Obtain R.
    M=1;
    if emd <= Emd
        gamdec_vec = linspace(min(i_XY),max(i_XY),100);
        for mm = 1:length(gamdec_vec)
            gamdec = gamdec_vec(mm);
            % M_test = 1 + 2*((Eie - emd)/(1-emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec)); % This is an old version including the intersection that gave a factor 1-emd
            M_test = 1 + 2*((Eie - emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
            if M_test > M
                M = M_test;
            end
        end
    end
    R_test = log2(floor(M))/n;
    if R_test > R
        R = R_test;
    end
end

end

function [R, Px_1,P1] = Achievability_Th1_pow_skew(n, rho, Efa, Emd, Eie)
%
% Achievability_Th1_pow_skew(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 2 for
% detection and decoding for the bi-AWGN channel. This bound skew the
% inputs and perform power allocation. Power allocation do not provide
% significant gains in the performed experiments while complicating the
% computation. Hence, this version is not used in the paper. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR [dB]
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% Px_1 = Prob[X=1]

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

% Definitions according to the paper:
i_XY_fcn = @(Px,rho, x, y) log(exp(sqrt(rho)*y.*x) ./ (Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) ;
r_XY_fcn = @(Px,rho, y) -log((Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) + rho/2 ;
% r_XY_fcn2 = @(Px, y) -log((Px * exp(y*sqrt(rho)+rho) + (1-Px) * exp(-y*sqrt(rho)-rho))) + rho/2 ;

% STEP 1: Skew input: try to find Px such that emd and efa are satisfied.
steplength = 0.01;
step_P1 = rho/2;
P1 = rho -step_P1;
R = 0; 
FLAG_rate= 0; 
while FLAG_rate == 0
    emd=1;%initital values
    P1 = P1 + step_P1;
    Px_1 = 0.5-steplength;
    if P1 > n*rho
        return;
    end
    P2 = (n*rho-P1)/(n-1);
    while emd > Emd
        r_Y_noise = 0;
        %     r_Y = 0;
        Px_1 = Px_1 + steplength;
        if Px_1 == 0.5
            %generate random variable r(Y). (Definition in the paper)
            Y = gaussian_rv(0,1,1,N);
            r_Y_noise = r_Y_noise + r_XY_fcn(Px_1,P1, Y);
        else
            %generate random variable r(Y). (Definition in the paper)
            for i = 1:n
                Y = gaussian_rv(0,1,1,N);
                if i==1
                    r_Y_noise = r_Y_noise + r_XY_fcn(Px_1,P1, Y);
                else
                    r_Y_noise = r_Y_noise + r_XY_fcn(Px_1,P2, Y);
                end
            end
        end
        %Obtain ccdf of the r.v.'s
        [cdf_R, Rx] = ecdf(-r_Y_noise); %use built in function to obtain empirical cdf
        Rx(1) = Rx(1) + 1e-8;
        gamma_fa = interp1(cdf_R, Rx, 1-Efa); % this is the minimum value of the threshold we can use
        
        %[cdf_R2, Rx2] = ecdf(-r_Y); %use built in function to obtain empirical cdf
        %Rx2(1) = Rx2(1) + 1e-8;
        %emd2 = interp1(Rx2,cdf_R2, gamma_fa); % This version requires many more samples
        emd = mean(exp(-r_Y_noise) .* (-r_Y_noise < gamma_fa)); % Includes a change of measure.
        
        if Px_1 >= 1
            R=0; break;
        end
    end
    
    % STEP 2: Find the maximum rate achievabable. Extra skewness of the input
    % allowed if it leads to an improved rate. (This extra step seems not
    % to help here)
%     R_temp = 0;
%     R_test = 0;
%     Px_1 = Px_1 - steplength;
%     while R_temp <= R_test
%         r_Y_noise = 0;
%         Px_1 = Px_1 + steplength;
%         if Px_1 > 1
%             Px_1 = 1;
%             break;
%         end
%         %generate random variable r(Y): (Definition in the paper)
%         if Px_1==0.5
%             Y = gaussian_rv(0,1,1,N);
%             r_Y_noise = r_Y_noise + r_XY_fcn(Px_1,P1, Y);
%         else
%             for i = 1:n
%                 Y = gaussian_rv(0,1,1,N);
%                 if i==1
%                     r_Y_noise = r_Y_noise + r_XY_fcn(Px_1,P1, Y);
%                 else
%                     r_Y_noise = r_Y_noise + r_XY_fcn(Px_1,P2, Y);
%                 end
%             end
%         end
%         [cdf_R, Rx] = ecdf(r_Y_noise); %use built in function to obtain empirical cdf
%         %Rx(1) = Rx(1) + 1e-8;
%         gamma_fa = interp1(cdf_R, Rx, Efa); % this is the minimum value of the threshold we can use
%         
%         emd = mean(exp(-r_Y_noise) .* (r_Y_noise > gamma_fa));
        
        % Generate information density i(X;Y): (Definition in the paper)
        i_XY = 0;
        for i = 1:n
            if i==1
                X = (-1).^(rand(N,1) > Px_1); %Pr[X=1]=Px
                Y = sqrt(P1)*X + gaussian_rv(0,1,1,N);
                i_XY = i_XY + i_XY_fcn(Px_1,P1, X, Y);
            else
                X = (-1).^(rand(N,1) > Px_1); %Pr[X=1]=Px
                Y = sqrt(P2)*X + gaussian_rv(0,1,1,N);
                i_XY = i_XY + i_XY_fcn(Px_1,P2, X, Y);
            end
        end
        % STEP 3: Obtain R.
        M=1;
        if emd <= Emd
            gamdec_vec = linspace(min(i_XY),max(i_XY),100);
            for mm = 1:length(gamdec_vec)
                gamdec = gamdec_vec(mm);
                M_test = 1 + 2*((Eie - emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
                if M_test > M
                    M = M_test;
                end
            end
        end
        R_temp = log2(floor(M))/n;

    if R_temp < R
        FLAG_rate = 1;
    else
        R = R_temp;
    end
end
end

function [R,P1] = Achievability_Th1_pow(n, rho, Efa, Emd, Eie)
%
% Achievability_Th1_pow(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 1 for
% detection and decoding for the bi-AWGN channel. This version do not skew
% the input but only perform power allocation to facilitate detection. This
% strategy seems to be much worse than skewing the input, so it is not
% included in the paper. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR [dB]
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% Px_1 = Prob[X=1]

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

% Definitions according to the paper:
i_XY_fcn = @(Px,rho, x, y) log(exp(sqrt(rho)*y.*x) ./ (Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) ;
r_XY_fcn = @(Px,rho, y) -log((Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) + rho/2 ;
% r_XY_fcn2 = @(Px, y) -log((Px * exp(y*sqrt(rho)+rho) + (1-Px) * exp(-y*sqrt(rho)-rho))) + rho/2 ;

% STEP 1: Power allocation for symmetric input: try to find P1 such that emd and efa are satisfied.
step_P1 = rho/2;
P1 = rho -step_P1;
Px_1 = 0.5;
R=0;
emd=1;%initital values
    
while emd > Emd
    P1 = P1 + step_P1;
    if P1 > n*rho
        return;
    end
    P2 = (n*rho-P1)/(n-1);

    
    %generate random variable r(Y). (Definition in the paper)
    Y = gaussian_rv(0,1,1,N);
    r_Y_noise = r_XY_fcn(Px_1,P1, Y);
    
    %Obtain ccdf of the r.v.'s
    [cdf_R, Rx] = ecdf(-r_Y_noise); %use built in function to obtain empirical cdf
    Rx(1) = Rx(1) + 1e-8;
    gamma_fa = interp1(cdf_R, Rx, 1-Efa); % this is the minimum value of the threshold we can use
    
    % -------------------------------------
    % This version requires many more samples:
    %[cdf_R2, Rx2] = ecdf(-r_Y); %use built in function to obtain empirical cdf
    %Rx2(1) = Rx2(1) + 1e-8;
    %emd2 = interp1(Rx2,cdf_R2, gamma_fa);
    %--------------------------------------
    emd = mean(exp(-r_Y_noise) .* (-r_Y_noise < gamma_fa)); % Includes a change of measure.
        
end
        
% Generate information density i(X;Y): (Definition in the paper)
i_XY = 0;
for i = 1:n
    if i==1
        X = (-1).^(rand(N,1) > Px_1); %Pr[X=1]=Px
        Y = sqrt(P1)*X + gaussian_rv(0,1,1,N);
        i_XY = i_XY + i_XY_fcn(Px_1,P1, X, Y);
    else
        X = (-1).^(rand(N,1) > Px_1); %Pr[X=1]=Px
        Y = sqrt(P2)*X + gaussian_rv(0,1,1,N);
        i_XY = i_XY + i_XY_fcn(Px_1,P2, X, Y);
    end
end
% STEP 3: Obtain R.
M=1;
if emd <= Emd
    gamdec_vec = linspace(min(i_XY),max(i_XY),100);
    for mm = 1:length(gamdec_vec)
        gamdec = gamdec_vec(mm);
        M_test = 1 + 2*((Eie - emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
        if M_test > M
            M = M_test;
        end
    end
end
R = log2(floor(M))/n;

end

function R = Achievability_X(n, rho, Efa, Emd, Eie)
%
% Achievability_X(n, rho, Efa, Emd, Eie):
% Computation of an achievability bound not used finally in the paper since 
% it was shown to be loose compared to the presented. For 
% detection and decoding for the bi-AWGN channel. We use QY=PY. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR (not in dB)
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end
    % STEP 1: Generate random variables (in a for loop to save memory)
    batches = 100;
    Ntmp = N/batches;
    i_XY = nan(1,N);
    for i=1:batches
        U = gaussian_rv(rho, rho, n, Ntmp);
        i_XY((i-1)*Ntmp + 1 : i*Ntmp) = n*log(2) - sum(log(1+exp(-2*U)), 2);
    end
    
    % To use the saddlepoint approximation to compute the FA
    Z = gaussian_rv(0, rho, 1, N);
    F = log(2) - log(1+exp(-2*Z));
    mgf = @(X, tau) mean(exp(X * tau));
    cgf_p = @(X, tau)  mean(X .* exp(tau * X)) / mgf(X,tau);
    cgf_pp = @(X,tau) (mean(X.^2 .* exp(tau*X))*mgf(X,tau) -  mean(X .* exp(tau * X))^2) / mgf(X,tau)^2;
    
    
    [cdf_i, ix] = ecdf(i_XY); %use built in function to obtain empirical cdf
    ix(1) = ix(1) + 1e-8; %make first element unique
    
    % STEP 2: Obtain the largest M s.t. all the constraints are satisfied
    alpha_d_vec = linspace(0,Emd, 50); %alpha is chosen so that emd is satisfied
    M = 1;
    for k = 1:length(alpha_d_vec)
        alpha_d = alpha_d_vec(k);
        gamma = interp1(cdf_i,ix,alpha_d);
        
        % The FA probability term is obtained using the saddlepoint approx.
        % find tau for saddlepoint approximation:
        tau = linspace(0,5,100);
        cgfp_tmp = nan(size(tau));
        for tt = 1:length(tau)
            cgfp_tmp(tt) = cgf_p(F, tau(tt));
        end
        tau = interp1(n*cgfp_tmp - gamma, tau, 0);
        
        cgf = log(mgf(F, tau));
        cgfp = cgf_p(F,tau);
        cgfpp = cgf_pp(F,tau);
        
        M_test = M;
        eie=0; efa=0;
        while eie <= Eie && efa <= Efa
            M_test = 2*M_test;
            eie = (M_test-1)/2*mean(exp(-i_XY).*(i_XY >= gamma))+alpha_d;%get eie
            efa = M_test * exp(n*(cgf - tau*cgfp)) * exp((n*cgfpp * tau^2) / 2) * qfunc(tau*sqrt(n*cgfpp));%get efa
        end
        if M_test > M
            M = M_test/2;
        end    
    end
   R = log2(M)/n;

end

function [R, Px_1] = converse_Th2(n, rho, Efa, Emd, Eie)
%
% converse_Th3(n, rho, Efa, Emd, Eie):
% Computation of the converse bound given in Theorem 2 for joint
% detection and decoding for the bi-AWGN channel. We use QY=PY. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR (not in dB)
% efa   = False alarm probability constraint
% emd   = Misdetection probability constraint
% eie   = Inclusive error probability constraint (includes the MD)
% OUTPUTS:
% R    = Rate (bits)

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

% Definitions according to the paper:
i_XY_fcn = @(Px, x, y) log(exp(y.*x) ./ (Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) ;
r_XY_fcn = @(Px, y) -log((Px * exp(y*sqrt(rho)) + (1-Px) * exp(-y*sqrt(rho)))) + rho/2 ;

% STEP 1: Skew input: Find Px such that emd and efa are satisfied.
steplength = 0.01;
Px_1 = 0.5-steplength;
emd=1;%initital values
while emd > Emd
    r_Y_noise = 0;
    Px_1 = Px_1 + steplength;
    %generate random variable r(Y). (Definition in the paper)
    for i = 1:n
        Y = gaussian_rv(0,1,1,N);
        r_Y_noise = r_Y_noise + r_XY_fcn(Px_1, Y);
        %r_Y = r_Y + r_XY_fcn2(Px_1, Y);
    end
    
    %Obtain ccdf of the r.v.'s
    [cdf_R, Rx] = ecdf(r_Y_noise); %use built in function to obtain empirical cdf
    Rx(1) = Rx(1) + 1e-8;
    
    gamma_fa = interp1(cdf_R, Rx, Efa); % this is the minimum value of the threshold we can use
    
    emd = mean(exp(-r_Y_noise) .* (r_Y_noise > gamma_fa)); % Using the change of measure trick
    
    if Px_1 >= 1
        R=0; return;
    end
end

% STEP 2: Find the maximum rate achievabable. Extra skewness of the input
% allowed.
R = 0;
R_test = 0;
Px_1 = Px_1 - steplength;
while R <= R_test
    r_Y_noise = 0;
    i_XY = 0;
    Px_1 = Px_1 + steplength;
    if Px_1 > 1
        Px_1 = 1;
        return;
    end
    %generate random variable r(Y): (Definition in the paper)
    for i = 1:n
        %X = sqrt(rho)*(-1).^(rand(1,N) > Px); %Pr[X=1]=Px
        Y = gaussian_rv(0,1,1,N);
        r_Y_noise = r_Y_noise + r_XY_fcn(Px_1, Y);
    end
    
    [cdf_R, Rx] = ecdf(r_Y_noise); %use built in function to obtain empirical cdf
    %Rx(1) = Rx(1) + 1e-8;
    gamma_fa = interp1(cdf_R, Rx, Efa); % this is the minimum value of the threshold we can use
    
    beta_fa = mean(exp(-r_Y_noise) .* (r_Y_noise > gamma_fa)); % Using the change of measure trick
    
    % STEP 3: Computation of M. We only get here if the previous condition
    % is not true.
    % Generate information density i(X;Y): (Definition in the paper)
    for i = 1:n
        X = sqrt(rho)*(-1).^(rand(N,1) > Px_1); %Pr[X=1]=Px
        Y = X + gaussian_rv(0,1,1,N);
        i_XY = i_XY + i_XY_fcn(Px_1, X, Y);
    end
    %Obtain ccdf of the r.v.'s
    [cdf_i, ix] = ecdf(i_XY);
    
    %obtain the threshold for the converse
    xi = interp1(cdf_i ,ix, Eie);
    
    %compute the converse
    beta_ie = mean(exp(-i_XY).*(i_XY >= xi));
    R_test = log2((1-beta_fa) / beta_ie) / n;
    %R_test = log2((1) / beta_ie) / n;
    if R_test > R
        R = R_test;
    end
end

    
end


function [R, npre] = Achievability_Preamble(n, rho, Efa, Emd, Eie)
%
% Achievability_Preamble(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 2 for
% detection and decoding for the bi-AWGN channel. We use QY=PY. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR (not in dB)
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% npre = optimal size of the preamble

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

for npre = 1:n-1
    thres = qfuncinv(Efa)*sqrt(npre * rho); %choose threshold to satisfy Efa
    emd = 1-qfunc((thres - npre*rho)/sqrt(npre*rho)); %what is the corresponding emd
    if emd <= Emd
        break;
    end
end

if emd > Emd
    R=0;
    return;
end

R = zeros(1,n-npre);
npre_pos = npre:n-1;
for nd = n-npre:-1:1
    npre = n - nd;
    thres = qfuncinv(Efa)*sqrt(npre * rho); %choose threshold to satisfy Efa
    emd = 1-qfunc((thres - npre*rho)/sqrt(npre*rho)); %what is the corresponding emd
    
    batches = 100;
    Ntmp = N/batches;
    i_XY = nan(1,N);
    for i=1:batches
        U = gaussian_rv(rho,rho,nd,Ntmp);
        i_XY((i-1)*Ntmp + 1 : i*Ntmp) = nd*log(2) - sum(log(1+exp(-2*U)),2);
    end
    
    M=1;
    if emd <= Emd
        gamdec_vec = linspace(min(i_XY),max(i_XY),100);
        for mm = 1:length(gamdec_vec)
            gamdec = gamdec_vec(mm);
            M_test = 1 + 2*((Eie - emd)/(1-emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
            %M_test = 1 + 2*((Eie - emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
            if M_test > M
                M = M_test;
            end
        end
    end
    R(nd) = log2(floor(M))/n;
end
R = fliplr(R);
[R,idx] = max(R);
npre = npre_pos(idx);
end

function [R,P1] = Achievability_Preamble_pow(n, rho, Efa, Emd, Eie)
%
% Achievability_Preamble(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 2 for
% detection and decoding for the bi-AWGN channel. We use QY=PY. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR (not in dB)
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% npre = optimal size of the preamble

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

P1_vec = linspace(rho,n*rho,300);


for pp = 1:length(P1_vec)
    P1 = P1_vec(pp);
    thres = qfuncinv(Efa)*sqrt(P1); %choose threshold to satisfy Efa
    emd = 1-qfunc((thres - P1)/sqrt(P1)); %what is the corresponding emd
    if emd <= Emd
        P2 = (n*rho-P1)/(n-1);
        break;
    end
end

if emd > Emd
    R=0;
    return;
end

nd = n -1;


batches = 100;
Ntmp = N/batches;
i_XY = nan(1,N);
for i=1:batches
    U = gaussian_rv(P2,P2,nd,Ntmp);
    i_XY((i-1)*Ntmp + 1 : i*Ntmp) = nd*log(2) - sum(log(1+exp(-2*U)),2);
end

M=1;
if emd <= Emd
    gamdec_vec = linspace(min(i_XY),max(i_XY),100);
    for mm = 1:length(gamdec_vec)
        gamdec = gamdec_vec(mm);
        M_test = 1 + 2*((Eie - emd)/(1-emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
        %M_test = 1 + 2*((Eie - emd) - mean(i_XY<gamdec))./mean(exp(-i_XY).*(i_XY>=gamdec));
        if M_test > M
            M = M_test;
        end
    end
end
R = log2(floor(M))/n;


end

function [R, npre] = Converse_Preamble(n, rho, Efa, Emd, Eie)
%
% Converse_Preamble(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 2 for
% detection and decoding for the bi-AWGN channel. We use QY=PY. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR (not in dB)
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% npre = optimal size of the preamble

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

for npre = 1:n-1
    gamma = sqrt(npre*rho)*qfuncinv((1-Emd)) + npre*rho; %find threshold based on Emd
    if qfunc(gamma/sqrt(npre*rho)) <= Efa %do we satisfy the false alarm with the threshold?
        break
    end
end

nd = n - npre; %remaining symbols used for data

if qfunc(gamma/sqrt(npre*rho)) > Efa
    R=0;
else
    iXY = 0;
    for i=1:nd
        U = gaussian_rv(rho, rho, 1, N);
        iXY = iXY + log(2) - log(1+exp(-2*U));
    end
    iXY = sort(iXY);
    %Obtain ccdf of the r.v.'s
    [cdf_i, ix] = ecdf(iXY);
    ccdf_i = 1-cdf_i;
    % Obtain the thresholds for the converse
    xi = interp1(ccdf_i ,ix, 1-Eie);
    % Obtain the rate
    R = log2(1 / mean(exp(-iXY).*(iXY >= xi))) / n;
end
end

function [R, P1] = Converse_Preamble_pow(n, rho, Efa, Emd, Eie)
%
% Converse_Preamble(n, rho, Efa, Emd, Eie):
% Computation of the achievability bound given in Theorem 2 for
% detection and decoding for the bi-AWGN channel. We use QY=PY. 
%
% INPUTS:
% n     = Blocklength
% rho   = SNR (not in dB)
% Efa   = False alarm probability constraint
% Emd   = Misdetection probability constraint
% Eie   = Inclusive error probability constraint (includes the MD)
%
% OUTPUTS:
% R    = Rate (bits)
% npre = optimal size of the preamble

N = round(100/min(min(Eie,Emd), Efa)); %number of samples of random variables
if N > 1e7
    N=1e7;
end

P1_vec = linspace(rho,n*rho,300);


for pp = 1:length(P1_vec)
    P1 = P1_vec(pp);
    thres = qfuncinv(Efa)*sqrt(P1); %choose threshold to satisfy Efa
    emd = 1-qfunc((thres - P1)/sqrt(P1)); %what is the corresponding emd
    if emd <= Emd
        P2 = (n*rho-P1)/(n-1);
        break;
    end
end

if emd > Emd
    R=0;
    return;
end

nd = n - 1; %remaining symbols used for data

iXY = 0;
for i=1:nd
    U = gaussian_rv(P2, P2, 1, N);
    iXY = iXY + log(2) - log(1+exp(-2*U));
end
iXY = sort(iXY);
%Obtain ccdf of the r.v.'s
[cdf_i, ix] = ecdf(iXY);
ccdf_i = 1-cdf_i;
% Obtain the thresholds for the converse
xi = interp1(ccdf_i ,ix, 1-Eie);
% Obtain the rate
R = log2(1 / mean(exp(-iXY).*(iXY >= xi))) / n;
end


function Z = gaussian_rv(mu,var, n, samples)

Z = sqrt(var) * randn(samples,n) + mu;

end

function R_bits = metaconverse_biawgn_fixed_s(eps,n,snr,s)
%
% R_bits = metaconverse_biawgn_fixed_s(eps_vec,n,snr,s,N):
% Computation of the metaconverse bound for the binary-input AWGN channel.
%
% INPUTS:
% eps = Probability of error
% n   = Blocklength
% snr = Value of the snr in linear scale (no dB!)
% s   = Parameter s of the generalized information density (s>0)
%
% OUTPUT:
% R_bits = Rate (bits) obtained as a result of the computation of the bound

N = round(100/eps); %number of samples of random variables
if N > 1e7
    N=1e7;
end


i_s=info_dens_biawgn(snr,n,s,N);
mu_f = @(y) (1/(2*sqrt(2*pi)^s) * (exp(-s/2*(y+sqrt(snr)).^2) + exp(-s/2*(y-sqrt(snr)).^2))).^(1/s);
Zmin = -10; Zmax = 10;
logmu = log(integral(mu_f, Zmin, Zmax));
j_s = n*logmu + i_s/s;
[cdf,gamma] = ecdf(j_s); % Compute the empirical cdf and save the values where it was evaluated (gamma).

mu = interp1(cdf ,gamma, eps);
R_bits = -1/n*log2(mean(exp(-j_s).*(j_s>= mu)));

end

function R = perfect_achievability_DT(n, snr, epsil)
%
% R = perfect_achievability_DT(n, snr, epsil):
% Computation of the DT bound.
%
% INPUTS:
% n     = Blocklength
% snr   = Value of the snr in linear scale (no dB!)
% epsil = probability of error
%
% OUTPUT:
% R = Rate (bits)


N = round(1000/epsil); %number of samples of random variables
if N > 1e7
    N=1e7;
end
TOL = epsil/50;
est_err = 1;
epsil_cur = 0;
R=0;
steplength = 0.3;
while est_err > TOL
    if epsil_cur > epsil
        R = R-steplength;
        steplength = steplength/2;
        R = R+steplength;
    else
        R = R+steplength;
    end
    
    if steplength < 2^-20
        R = 0; break;
    end
    
    i_s = info_dens_biawgn(snr,n,1,N);
    epsil_cur = mean( exp( - max(0, i_s - log((exp(n*R*log(2))-1)/2))));
    
    est_err = abs(epsil_cur - epsil);
end
end

function C = capacity(rho)

sigma_sq = 1;
Omega = rho/sigma_sq;
h = @(x) log(1+exp(-2*x));
intgrl_fcn = @(x,l, s) exp(-(x-Omega).^2 / (2*Omega)).*exp(-s*h(x)).*(-h(x)).^l;
H = @(l, s) (1/sqrt(2*pi*Omega)) * integral(@(x) intgrl_fcn(x,l,s), -100,100);

H1 = H(1,0);
C = 1+H1/log(2); %use Erseghe's method to compute capacity

end

function [n, out_val] = NA(rho,n,bler_target)
% This approximation is based on eq. 95 in "?Coding in the Finite-Blocklength Regime: Bounds Based on Laplace Integrals and Their
%Asymptotic Approximations" By Tomaso Erseghe.


sigma_sq = 1;
Omega = rho/sigma_sq;

h = @(x) log(1+exp(-2*x));
intgrl_fcn = @(x,l, s) exp(-(x-Omega).^2 / (2*Omega)).*exp(-s*h(x)).*(-h(x)).^l;
H = @(l, s) (1/sqrt(2*pi*Omega)) * integral(@(x) intgrl_fcn(x,l,s), -100,100);

H1 = H(1,0);
H2 = H(2,0);

C = 1+H1/log(2);
V = H2-H1^2;
alpha = 1-1./sqrt(n);

out_val = max(C - sqrt(V./n)*log2(exp(1)).*qfuncinv(alpha*bler_target) + log2(1-alpha)./n + log2(n)./(2*n),0);



end

