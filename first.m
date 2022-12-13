[sl,sr] = gen_echo(1000, 30, 0, 10, 20000, 80000, 4);
[sf, slh, srh] = crt_sf(sl, sr, 80000, 4);

[sl2,sr2] = gen_echo(100, 45, 10, 5, 20000, 80000, 4);
[sf2, slh2, srh2] = crt_sf(sl2, sr2, 80000, 4);

base = 0.03;
c = 1500;
fs = 20e3;
index = [109801 111271 111904 112045 112537 113031 113666 114373 114515 115649 116076];
index2 = [13456 14877 15502 15642 16135 16632 17277 18000 18145 19317 19760];
[dfi, doa] = crt_dfi(slh, srh, index, base, fs, c);
[dfi2, doa2] = crt_dfi(slh2, srh2, index2, base, fs, c);

figure()
scatter(dfi,doa,'filled') 
grid on;
xlabel('dfi')
ylabel('doa')
hold on
scatter(dfi2,doa2,'filled') 
legend("1 цель", '2 цель')

function [dfi, doa] = crt_dfi(slh, srh, index, base, fs, c)
dfi=angle(slh(index))-angle(srh(index));
dt=dfi/(2*pi*fs);
doa=asin(c*dt/base);
end

function [sf, slh, srh] = crt_sf(sl,sr,fd,ts);
t = 0:1/fd:ts;
s = hilbert(gen_M_signal(t));
slh = hilbert(sl);
srh = hilbert(sr);
sf = xcorr(slh, s);
sf = sf(length(slh):end);
treshhold = zeros(size(sf));
for i = 1:100:(length(sf)-fd);
    skvo = sqrt(sum(abs(sf(i:i+fd)).^2)/fd);
    treshhold(fd/2+i:fd/2+i+99) = 4.8*skvo;
end
treshhold(1:40000) = 4.8*skvo;
treshhold(i:end) = 4.8*skvo;

figure()
plot(abs(sf));
hold on;
plot(treshhold, 'r');
end

function s = gen_M_signal(t);
F = 20e3;
Td = 0.001;
m = [0 0 0 0 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1 0 1 0 0 0 1 1 1 ...
    1 0 0 1 0 0 0 1 0 1 1 0 0 1 1 1 0 1 0 1 0 0 1 1 1 1 1 ...
    0 1 0 0 0 0 1 1 1 0 0 0 1 0 0 1 0 0 1 1 0 1 1 0 1 0 1 ...
    1 0 1 1 1 1 0 1 1 0 0 0 1 1 0 1 0 0 1 0 1 1 1 0 1 1 1 ...
    0 0 1 1 0 0 1 0 1 0 1 0 1 1 1 1 1 1 1];
len = length(t);
s = zeros(1,len);
if length(t) == 1
    if (t>0)&(t<0.127)
        sign = m(floor(t*1000)+1)*2-1;
        s = sign*sin(2*pi*F*t);
    else
        s = 0;
    end
else
    for i = 1:length(t)
        if (t(i)>0)&(t(i)<0.127)
            sign = m(floor(t(i)*1000)+1)*2-1;
            s(i) = sign*sin(2*pi*F*t(i));
        end
    end
end
end

function [sl, sr] = gen_echo(R_distance, R_phi, KUt, Vt, fs, fd, T)
num = floor(T*fd)+1;
c = 1500;
base = 0.03;
lens = [100 110 91 84 75 66 77 112 45 134 128]/1000*c/2;
amps = [1.3 1.8 2.9 7.5 6.8 1.3 4.1 3.2 2.3 1.8];
sl = zeros(1,num);
sr = zeros(1,num);
t = 0:1/fd:T;
for i=1:length(lens)
    xt = R_distance*cos(R_phi/180*pi) + lens(i)*cos(KUt/180*pi);
    yt = R_distance*sin(R_phi/180*pi) + lens(i)*sin(KUt/180*pi);
    R_dist_t = sqrt(xt.^2+yt.^2);
    R_phi_t = atan(yt/xt);
    fs_t = fs*(1-2*Vt/c);
    sl = sl + gen_M_signal(t-2*R_dist_t/c);
    sr = sr + gen_M_signal(t-2*R_dist_t/c-base*sin(R_phi_t)/c);
end

common = randn(size(sl));
power = 1;
corr = 1/2;
sl = sl+power*(corr*common+(1-corr)*randn(size(sl)));
sr = sr+power*(corr*common+(1-corr)*randn(size(sl)));
end
