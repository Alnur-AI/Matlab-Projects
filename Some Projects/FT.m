%%

%% Задание 10
%funcs
func1 = @(t) sin(3*t)./(2*t);
func2 = @(t) 0+cos(t).*(t>=-1/8 & t<=1/8);
func3 = @(t) exp(-2*abs(t)).*log(1+t.^4);
func4 = @(t) t.^2.*exp(-t.^4);

func5 = @(L) 0+pi/2*(L<=3 & L>=-3);
func6 = @(L) 1./L;
func7 = @(t) (cos(t) - exp(-abs(t)))./t;

ftfunc1 = @(L) 0+pi/2*(L<=3 & L>=-3);
ftfunc2 = @(L) 0+1i*(sin((1+L)./3)./(1+L)-sin((1-L)./3)./(1-L));
%figure
figure1 = figure;

%figure2 = figure; % + оси
%steps
step0 = 0.1;
step1 = 0.1;
step2 = 0.05;
step3 = 0.01; 
step4 = 0.001;
step5 = 0.0001;
%inpLimVec
inpLimVec1 = [-3 10];
inpLimVec2 = [-100 100];
inpLimVec3 = [-4 4];
%outLimVec
outLimVec1 = [-10 10];
outLimVec2 = [-10 10];
outLimVec3 = [-30 30];
plotFT(figure1, func7, [], step5, [-100 90], outLimVec1);

function res = plotFT(hFigure, fHandle, fFTHandle, step, inpLimVec, outLimVec)
    SPlotInfo=get(hFigure,'UserData');
    if (~isempty(SPlotInfo))
        if (isempty(outLimVec))
            outLimVec = SPlotInfo.outLimVec;
        end
    end
    if (isempty(outLimVec))
        outLimVec = [-10 10];
    end
    figure(hFigure);
    N = round((inpLimVec(2) - inpLimVec(1))/step);
    step_new = (inpLimVec(2) - inpLimVec(1))/N;
    y = 0;
    if (isempty(fFTHandle)==0)
        x = linspace(outLimVec(1),outLimVec(2),1000);
        y = fFTHandle(x);
        hold on
        subplot(2,1,1);
        plot(x,real(y));
        subplot(2,1,2);
        plot(x,imag(y));
    end
    t=linspace(0.0000000001,inpLimVec(2)-inpLimVec(1),N);
    f = zeros(1,N);
    sm = inpLimVec(1)/(inpLimVec(2)-inpLimVec(1));
    c = ceil(sm);
    i = 1;
    while t(i)+c*(inpLimVec(2)-inpLimVec(1))<=inpLimVec(2)
        f(i) =  fHandle(t(i)+c*(inpLimVec(2)-inpLimVec(1)));
        i =i+1;
    end
    add = step_new - (-inpLimVec(2)+t(i)+c*(inpLimVec(2)-inpLimVec(1)));
 
    j = 0;
    while i+j<=N
        f(i+j) =  fHandle(t(j+1)-step_new+inpLimVec(1)+add);%!
        j =j+1;
    end
    T = abs(inpLimVec(2) - inpLimVec(1));
    F = T*fft(f)/(N);
    Nyq=2*pi*N/(T);
    df=2*pi/T;
    nu=-5*Nyq+df*(0:10*N);
    F_2 = [F F F F F F F F F F];
    subplot(2,1,1);
    axis ([outLimVec(1) outLimVec(2) min(min(real(F)),min(real(y))) max(max(real(F)),max(real(y)))]);
    hold on 
    pRE = plot(nu(1:10*N), (real(F_2(1:10*N))),'r');
    legend('FT','approx FT');
    xlabel('L');
    ylabel('Re(L)');
    nu=-5*Nyq+df*(0:10*N);
    subplot(2,1,2);
    hold on 
    xlabel('L');
    ylabel('Imag(L)');
    axis ([outLimVec(1) outLimVec(2) min(min(imag(F)),min(imag(y))) max(max(imag(F)),max(imag(y)))]);
    pIM = plot(nu(1:10*N), (imag(F_2(1:10*N))),'r');
    legend('FT','approx FT');
    str = struct('Re_hangle',@pRE,'Im_hangle',@pIM,'outLimVec',outLimVec); 
    if (isempty(SPlotInfo))
        SPlotInfo = str;
    else
        SPlotInfo(end)= str;
    end
    set(hFigure, 'UserData',SPlotInfo);
    res = struct('nPoints',N, 'step', step_new ,'inpLimVec', inpLimVec,'outLimVec',outLimVec);
    disp(res);
end