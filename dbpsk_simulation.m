clear; clc; close all;

aa = 1;   %%runs
xx=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_symb = 5;                             %number of symbols sent
%         bit = randi([0,1],1,num_symb);
bit = [0 1 0 1 0];
LB = length(bit);                   %signal length
T = 10^-8;                          %bit duration
bw = 1/T;                           %bit rate and bw
f = 2/(2*T);                        %frequency
Fs = 3.02e+9;                       %sample frequency
Ts = 1/Fs;                          %sample rate
sps = ceil(T/Ts);                   %samples per second // baud rate
t_mod = 0:Ts:T;
beta = 0.5;                      %roll off factor
span = 10;                          %filter span in [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sampling Intervals%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
packet_duration = (0:Ts:LB*T-Ts);
bit_duration = (0:Ts:T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Signals%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sig1 = sin(2*pi*f*t_mod);
Sig2 = sin(2*pi*f*t_mod + pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase Noise %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;                      %duration
N = 465;                    %N samples
dt = t/(length(N));
dW = zeros(1,length(N));    %allocated space for dW
W = zeros(1,length(N));     %allocate space for W

%phase noise
phn = 10e-4; %phase noise variance
dW = sqrt(phn).*(randn(N,1));
W(1) = dW(1);
for j = 2:N
    dW = sqrt(phn).*(randn(N,1));           %real valued normal distribution AWGN
    W(j) = W(j-1) + dW(j);        %
end
v = var(dW); % 19e-4

EbN0dB = 1:5:31;
for ss = 1:5:31
    xx=xx+1;
    for ii=1:aa
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Modulation%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mod = [];                       %array for modulated signal
        pulse = [];

        for n = 1:LB
                pulse = zeros(1,sps);
                pulse(:,1) = bit(:,n)*2-1;
            
        mod = [mod pulse];
        end
        
%         mod = mod;
        t_mod = 0:Ts:num_symb*(ceil(T/Ts))*Ts-Ts;  %modulated signal timing

%         figure()
%         plot(t_mod,mod)               
%         title('Modulated Signal')
%         xlabel('Time [s]')
%         xlim([0 5.1e-8])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rasied Cosine Filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rc = rcosdesign(beta, span, sps, 'normal');  %raised cosine filter
%         rc_timing = (0:Ts:Ts*(length(rc)-1));

%         figure()
%         plot(rc)
%         xlim([0 311])
%         title("Raised Cosine Filter")
%         xlabel("Samples")
%         ylabel("Normalized Amplitude")
        
        pulse_shape = conv(mod.', rc.');
        ps = length(pulse_shape);
        t_shape = 0:Ts:(ps-1)*Ts;
%         figure()
%         plot(t_shape, pulse_shape)
%         title("Pulse Shaped Modulated Signal")
%         xlabel("Samples")
%         ylabel("Amplitude")
%         xlim([0 465])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Upconversion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         upconv = pulse_shape.' .* (cos(2*pi*(73e+9)*t_shape) + sin(2*pi*(73e+9)*t_shape)).';
        

%         figure()
%         plot(t_shape, upconv)
%         title("Upconverted Signal")
%         xlabel("Time [s]")
%         ylabel("Amplitude")tx = mod + noise.';
%         grid on

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Noise%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EbN0_linear = 10.^((ss)/10);
        Eb = sum(abs(pulse_shape).^2)/(length(pulse_shape));     %compute energy per bit
        N0 = Eb/EbN0_linear;                     %required NSD from Eb/N0
        noise = sqrt(N0)*(randn(length(pulse_shape),1)+1j*randn(length(pulse_shape),1));
        ber_TH = berawgn(EbN0dB, 'psk',2,'nondiff');
        
        
        new_sig = pulse_shape .* exp(1j.*W).';
        new = new_sig + noise;

        tx = pulse_shape.' + noise.';
        
     
%         figure()
%         plot(t_shape,tx)
%         title("Pulse Shaped Signal Without Phase Noise")
%         xlabel("Time [s]")
%         ylabel("Amplitude")      
%         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Fast Fourier Transform %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         D = fft(I);
        %         E = fftshift(D);
        %         f = (0:length(D)-1)*Fs/length(D);
        %         % plot(f,abs(E))
        %         fshift = (-Ld/2:Ld/2-1)*(Fs/Ld);
        %         powershift = abs(E).^2/Ld;

        % figure()
        % plot(fshift,abs(E))
        % xlabel('Frequency')
        % ylabel('Magnitude [dB]')
        % title('Modulated Signal Frequency Domain')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%tx
        %% Received Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         received = upconv + noise.';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Lowpass Filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         lp = lowpass(received, 100e+6, Fs)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Introduce Frequency/Phase Offset
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f_offset = 73e+9;
        p_offset = 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Downconversion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         received1 = s; 
        
%         figure()
%         plot(t_shape, new)
%         title("Received Signal with Phase Noise")
%         xlabel("Time [s]")
%         ylabel("Amplitude")
%         grid on

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Raised Cosine Filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rc0 = rcosdesign(beta, span, sps, 'normal');  %raised cosine filter
%        rc_timing = (0:Ts:Ts*(length(rc)-1));
        
%         figure()
%         plot(rc0)
%         xlim([0 311])
%         title("Raised Cosine Filter")
%         xlabel("Samples")
%         ylabel("Normalized Amplitude")
        
        pulse_shape0 = conv(new, rc0); % pulse shaping of signal with phase noise
        ps0 = length(pulse_shape0);
        t_shape0 = 0:Ts:(ps0-1)*Ts;
        
        pulse_shape1 = conv(tx, rc0); % pulse shaping of signal without phase noise
        ps1 = length(pulse_shape1);
        t_shape1 = 0:Ts:(ps1-1)*Ts;

%         figure()
%         plot(t_shape0, pulse_shape0)
%         title("Pulse Shaped Downconverted Signal")
%         xlabel("Time [s]")
%         ylabel("Amplitude")
%         xlim([0 2.6e-7])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Lowpass Filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lp_phn = lowpass(pulse_shape0, 100e+6, Fs); % filter modulated signal with phase noise
% 
%         figure()
%         plot(t_shape0,lp_r)
%         title("Receieved Filtered Signal Phase Noise")
%         xlabel("Time [s]")
%         ylabel("Amplitude")
%         grid on
        
        lp_awgn = lowpass(pulse_shape1, 100e+6, Fs); % filter mod signal with AWGN

%         figure()
%         plot(t_shape0,lp_r)
%         title("Receieved Filtered Signal AWGN")
%         xlabel("Time [s]")
%         ylabel("Amplitude")
%         grid on
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Matched Filters%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mark

        inst_sample = lp_phn(302:sps:end); % filter mod signal with phase noise
        inst_sample1 = lp_awgn(302:sps:end); % filter mod signal with AWGN

        % matched0_plot = 1:numel(matched0);

        % space
        % matched1 = conv(received, fliplr(Sig2));
        % inst_sample1 = abs(matched1(sps:sps:end-sps));
        % matched1_plot = 1:numel(matched1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Demodulator with Phase Noise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        demod0 = [];
        BER0 = [];

        yy = 1;
        for l = 300:sps:length(lp_phn)-350
        
            if lp_phn(l) < 0
                PP = 0;
            else
                PP = 1;
            end

            yy = yy+1;
            demod0 = [demod0 PP];
        end
        demod0

        err0 = 0;
        for c = 1:LB
            if bit(c) ~= demod0(c)
                err0 = err0+1;
            end
            BER0 = err0/length(bit);
        end
        BER1_0(ii) = BER0;
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Demodulator without Phase Noise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        demod1 = [];
        BER1 = [];

        XX=1;
        for k = 300:sps:length(lp_awgn)-350
        
            if lp_awgn(k) < 0
                OO = 0;
            else
                OO = 1;
            end
            XX = XX+1;
            demod1 = [demod1 OO];
        end
        demod1

        err1 = 0;
        for a = 1:LB
            if bit(a) ~= demod1(a)
                err1 = err1+1;
            end
            BER1 = err1/length(bit);
        end
        BER1_1(ii) = BER1;

    end     %iteration end
    AVG_BER0(xx)=sum(BER1_0)/(aa);
    AVG_BER1(xx)=sum(BER1_1)/(aa);
end         %end of Eb/N0

EbN0dB = 1:5:31;
figure()
semilogy(EbN0dB,AVG_BER0)
hold on
semilogy(EbN0dB,AVG_BER1)
grid on
legend("BER, Phase Noise", "BER, AWGN")
title("Comparison of BER with Phase Noise and without")
xlabel('Eb/No (dB)')
ylabel('BER')
