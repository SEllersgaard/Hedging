%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This programme computes and plots the optimal market surfaces 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global K;
global vol;
global T;
global N;
global Y;
global kappa;
global lambda;
global g;
global Jt;
global dt;
global ds;
global smin;
global smax;
global deltaPLUS;
global deltaMINUS;
global Upper;
global Lower;
global Mid;
global probability;
global cut;
global qmin;


K   = 50.5;     % Strike
vol = 0.5;      % Volatility
T   = 1;        % Initial Time To Maturity
N   = 100;      % Number of calls sold
Y   = 1;        % Surcharge From Market Order (vis-a-vis half-spread)
Y2  = Y;
phi = 1;      % Risk Aversion to Being Off Perfect Hedge Surface
kappa = 0.3;    % Decay Factor in Arrival Orders 
lambda = 100;  %50;    % Poisson Intensity
g      = 200;    % Poisson Intensity Two

%%%%%%%%% Allow for jumps to non-adjacent nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%
probability   = 0.90;
cut = 4;
prob=zeros(1,cut);
for i=1:1:cut-1
   prob(i) = (1-probability)^(i-1)*probability; 
end
prob(cut) = 1-sum(prob(1:cut-1));

vol2   = sqrt(2*g)*(sqrt(2-probability)/probability)*vol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


qmax = 120;
qmin = -20;
smax = 100;
smin = 0;

Js   = (smax-smin)*(1/vol);
Jt   = 1000; %1000
Jq   = qmax-qmin;
Jq2  = qmax;
ds   = (smax-smin)/Js;
dt   = T/Jt;
dq   = 1;

grid = zeros(Jq+1,Js+1,Jt+1);
marketgrid = zeros(Jq+1,Js+1,Jt);
optimalgrid = zeros(Jq+1,Js+1,Jt);


% Fill out the terminal condition

for s = 1:1:Js+1
   for q = 1:1:Jq+1 
    
       S = (smin+(s-1)*ds);
       Q = (qmin + (q-1)*dq);
       
       if Q >= N && S >= K
           grid(q,s,Jt+1) = N*K + (Q-N)*(S-Y2);
       elseif Q >= N && S < K
           grid(q,s,Jt+1) = Q*(S-Y2);
       elseif Q < N && S >= K
           grid(q,s,Jt+1) = N*(K) - (N-Q)*(S+Y);
       elseif Q < N && S < K
           grid(q,s,Jt+1) = Q*(S-Y2);
       end
       
   end
end

deltaMINUS = zeros(Jq+1,Js+1,Jt);
sellLimit = zeros(Jq+1,Js+1,Jt);
sellMarket = zeros(Jq+1,Js+1,Jt);

deltaPLUS = zeros(Jq+1,Js+1,Jt);
buyLimit = zeros(Jq+1,Js+1,Jt);
buyMarket = zeros(Jq+1,Js+1,Jt);

tic
% Go backwards in time (Explicitly)
for t = Jt:-1:1
    
    ttm = (Jt-(t-1))*dt;
    
    for s = 2:1:Js
        
        S = (smin+(s-1)*ds);
        d1 = (S-K)/(vol2*sqrt(ttm));
        Delta = N*normcdf(d1);
        %d1 = (log(S/K) + 0.5*vol^2*ttm)/(vol*sqrt(ttm));
        %Delta = N*normcdf(d1);
        
        DeltaRound = round(Delta,0);
        if DeltaRound == 0
            DeltaRound = 1;
        elseif DeltaRound == Jq + 2
            DeltaRound = Jq + 1;
        end
               
        
        for q = 2:1:Jq
        
            
            if s > cut && s <= Js+1-cut
                    sumLeft = 0;
                    sumRight = 0;
                    for l = 1:1:cut
                       sumLeft = sumLeft + prob(l)*grid(q,s-l,t+1);
                       sumRight = sumRight + prob(l)*grid(q,s+l,t+1);
                    end
            elseif s <= cut
                                    
                probL=zeros(1,s-1);
                if s > 2
                    for i=1:1:s-2
                       probL(i) = (1-probability)^(i-1)*probability; 
                    end
                    probL(s-1) = 1-sum(probL(1:s-2));
                else
                    probL(1) = 1;
                end
                sumLeft = 0;
                sumRight = 0;
                for l = 1:1:s-1
                   sumLeft = sumLeft + probL(l)*grid(q,s-l,t+1);
                end
                %sumLeft = grid(q,s-l,t+1);
                for l = 1:1:cut
                   sumRight = sumRight + prob(l)*grid(q,s+l,t+1);
                end
                
            elseif s > Js+1-cut
                s00 =  s - (Js+1-cut);
                probR=zeros(1,cut-s00);
                if s < Js
                    for i=1:1:cut-s00-1
                       probR(i) = (1-probability)^(i-1)*probability; 
                    end
                    probR(cut-s00) = 1-sum(probR(1:cut-s00-1));
                else
                    probR(1) = 1;
                end
                sumLeft = 0;
                sumRight = 0;
                for l = 1:1:cut
                   sumLeft = sumLeft +  prob(l)*grid(q,s-l,t+1);
                end
                for l = 1:1:cut-s00
                   sumRight = sumRight + probR(l)*grid(q,s+l,t+1);
                end
                %sumRight = grid(q,s+l,t+1);
            end
            
            
            
            
            Q = (qmin + (q-1)*dq);
           
            
            if Q >= Delta
                
                deltaMINUS(q,s,t) = 1/kappa - S + grid(q,s,t+1)-grid(q-1,s,t+1);
                    
                sellLimit(q,s,t) = grid(q,s,t+1) + dt*( ...
                        g*(sumLeft + sumRight - 2*grid(q,s,t+1)) + ...
                        (lambda/kappa)*exp(-kappa*deltaMINUS(q,s,t)) - ...
                        phi*(Q-Delta)^2);    
                                   
                sellMarket(q,s,t) = S - Y + grid(q-1,s,t+1);
                
                if sellMarket(q,s,t) >= sellLimit(q,s,t)
                   marketgrid(q,s,t) = -1; 
                end
                
                grid(q,s,t) = max(sellLimit(q,s,t),sellMarket(q,s,t));
                
            else
               
                deltaPLUS(q,s,t) = 1/kappa + S + grid(q,s,t+1)-grid(q+1,s,t+1);
                        
                buyLimit(q,s,t) = grid(q,s,t+1) + dt*( ...
                    g*(sumLeft + sumRight -2*grid(q,s,t+1)) + ...
                    (lambda/kappa)*exp(-kappa*deltaPLUS(q,s,t)) - ...
                    phi*(Q-Delta)^2);
                
                buyMarket(q,s,t) = -(S + Y) + grid(q+1,s,t+1);
                
                if buyMarket(q,s,t) >= buyLimit(q,s,t)
                   marketgrid(q,s,t) = 1; 
                end
                
                grid(q,s,t) = max(buyLimit(q,s,t),buyMarket(q,s,t));
                
            end
               
            
        end % end q loop
        
        marketgrid(-qmin+DeltaRound,s,t) = 2;
        optimalgrid(-qmin+DeltaRound,s,t) = 2;
        
    end % end s loop
    
    %Fill out ends of market grid
    marketgrid(Jq+1,:,t) = marketgrid(Jq,:,t);
    marketgrid(1,:,t) = marketgrid(2,:,t);
    marketgrid(:,Js+1,t) = marketgrid(:,Js,t);
    marketgrid(:,1,t) = marketgrid(:,2,t);
    
    %Extrapolate the ends of the grid
    grid(Jq+1,:,t) = grid(Jq,:,t) + (grid(Jq,:,t)-grid(Jq-1,:,t));
    grid(1,:,t) = grid(2,:,t) - (grid(3,:,t)-grid(2,:,t));
    grid(:,Js+1,t) = grid(:,Js,t) + (grid(:,Js,t)-grid(:,Js-1,t));
    grid(:,1,t) = grid(:,2,t) - (grid(:,3,t)-grid(:,2,t));
    
    
    
    
    
end % end t loop
toc

%marketgrid2 = zeros(qmax+1,Js+1,Jt);

marketgrid2 = marketgrid(-qmin+1:end,:,:);





Upper = zeros(Jt,Js+1);
Lower = zeros(Jt,Js+1);

for t = Jt:-1:1
   
    for s=1:1:Js+1
    
        for q=1:1:Jq2+1
           
            if marketgrid2(q,s,t) == -1
                Upper(t,s) = (q-1)*dq;
                if Upper(t,s) >= qmax
                    Upper(t,s) = qmax;
                end
                break;
            end
                
            if q == Jq2+1
                Upper(t,s) = qmax;
            end
            
        end
        
         for q=Jq2+1:-1:1
           
            if marketgrid2(q,s,t) == 1
                Lower(t,s) = (q-1)*dq;
                %if Lower(t,s) == -1
                %    Lower(t,s) = -1;
                %end
                break;
            end
                
            if q == 1
                Lower(t,s) = NaN;
            end
            
        end
      
    end
    
end

%Upper = flipud(Upper);
%Lower = flipud(Lower);

Mid = zeros(Jt,Js+1);
for t = 0:1:Jt-1;
   ttm = T-t*dt;
   for s = 0:1:Js
      
        S = (smin+(s-1)*ds);
        d1 = (S-K)/(vol2*sqrt(ttm));
        Mid(t+1,s+1) = N*normcdf(d1);
       
         %d1 = (log((s+1)*ds/K) + 0.5*vol^2*ttm)/(vol*sqrt(ttm));  
         %Mid(t+1,s+1) = N*normcdf(d1);
             
   end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin = 50;                 %Select how close to maturity we should go (0 for maturity)
time = T-dt:-dt:tmin*dt;
stock = smin:ds:smax;
low = 50;                   % Select lower range for stock dimension
high = 150;                 % Select upper range for stock dimension
topalpha = 1;
botalpha = 1;
topalpha2 = 1;             % Transparency of top surface
botalpha2 = 1;             % Transparency of bottom surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal SELL and MID surface

figure(1)

h1 = surf(stock(low:high),time,Upper(1:end-tmin,low:high));
set(h1,'LineStyle','none')
set(h1,'FaceAlpha',topalpha)
colormap(bone)
freezeColors
ylabel('TTM   ','FontSize',16)
xlabel('Stock price   ','FontSize',16)
zlabel('Inventory Level   ','FontSize',16)
zlim([0,100])
xlim([stock(low),stock(high)])
title('Sell Market Surface')

figure(2)
h2 = surf(stock(low:high),time,Lower(1:end-tmin,low:high));
set(h2,'LineStyle','none')
set(h2,'FaceAlpha',botalpha)
colormap(flipud(colormap(summer)))
freezeColors
ylabel('TTM   ','FontSize',16)
xlabel('Stock price   ','FontSize',16)
zlabel('Inventory Level  ','FontSize',16)
zlim([0,100])
xlim([stock(low),stock(high)])
title('Buy Market Surface')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(3)
subplot(2,2,1)

h1 = surf(stock(low:high),time,Upper(1:end-tmin,low:high));
set(h1,'LineStyle','none')
set(h1,'FaceAlpha',topalpha)
set(h1,'FaceColor',[53/255 42/255 134/255],'FaceAlpha',topalpha)
ylabel('TTM   ','FontSize',16)
xlabel('Stock price   ','FontSize',16)
zlabel('Inventory Level   ','FontSize',16)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(a)  ')

subplot(2,2,2)

h3 = surf(stock(low:high),time,Mid(1:end-tmin,low:high));
set(h3,'LineStyle','none')
set(h3,'FaceAlpha',topalpha)
colormap(autumn)
ylabel('TTM   ','FontSize',16)
xlabel('Stock price   ','FontSize',16)
zlabel('Inventory Level   ','FontSize',16)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(b)   ')

subplot(2,2,3)

h2 = surf(stock(low:high),time,Lower(1:end-tmin,low:high));
set(h2,'LineStyle','none')
set(h2,'FaceAlpha',botalpha)
%colormap(flipud(colormap(summer)))
set(h2,'FaceColor',[166/255 191/255 106/255],'FaceAlpha',topalpha2)
freezeColors
ylabel('TTM   ','FontSize',16)
xlabel('Stock price   ','FontSize',16)
zlabel('Inventory Level  ','FontSize',16)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(c)   ')

subplot(2,2,4)

h1 = surf(stock(low:high),time,Upper(1:end-tmin,low:high));
set(h1,'LineStyle','none')
set(h1,'FaceAlpha',topalpha2)
set(h1,'FaceColor',[53/255 42/255 134/255],'FaceAlpha',topalpha2)
hold on
h2 = surf(stock(low:high),time,Lower(1:end-tmin,low:high));
set(h2,'LineStyle','none')
set(h2,'FaceAlpha',botalpha2)
set(h2,'FaceColor',[166/255 191/255 106/255],'FaceAlpha',topalpha2)
hold on
h3 = surf(stock(low:high),time,Mid(1:end-tmin,low:high));
set(h3,'LineStyle','none')
set(h3,'FaceAlpha',topalpha2)
colormap(autumn)
ylabel('TTM   ','FontSize',16)
xlabel('Stock price   ','FontSize',16)
zlabel('Inventory Level   ','FontSize',16)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(d)  ')
%view(0,90)

figure(4)



subplot(1,3,1)

h1 = surf(stock(low:high),time,Upper(1:end-tmin,low:high));
set(h1,'LineStyle','none')
set(h1,'FaceAlpha',topalpha)
set(h1,'FaceColor',[53/255 42/255 134/255],'FaceAlpha',topalpha)
ylabel('TTM   ','FontSize',12)
xlabel('Stock Price   ','FontSize',12)
zlabel('Inventory Level   ','FontSize',12)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(a)   ')


subplot(1,3,2)

h2 = surf(stock(low:high),time,Lower(1:end-tmin,low:high));
set(h2,'LineStyle','none')
set(h2,'FaceAlpha',botalpha)
%colormap(flipud(colormap(summer)))
set(h2,'FaceColor',[166/255 191/255 106/255],'FaceAlpha',topalpha2)
freezeColors
ylabel('TTM   ','FontSize',12)
xlabel('Stock Price   ','FontSize',12)
zlabel('Inventory Level  ','FontSize',12)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(b)   ')

subplot(1,3,3)

h1 = surf(stock(low:high),time,Upper(1:end-tmin,low:high));
set(h1,'LineStyle','none')
set(h1,'FaceAlpha',topalpha2)
set(h1,'FaceColor',[53/255 42/255 134/255],'FaceAlpha',topalpha2)
hold on
h2 = surf(stock(low:high),time,Lower(1:end-tmin,low:high));
set(h2,'LineStyle','none')
set(h2,'FaceAlpha',botalpha2)
set(h2,'FaceColor',[166/255 191/255 106/255],'FaceAlpha',topalpha2)
hold on
h3 = surf(stock(low:high),time,Mid(1:end-tmin,low:high));
set(h3,'LineStyle','none')
set(h3,'FaceAlpha',topalpha2)
colormap(autumn)
ylabel('TTM   ','FontSize',12)
xlabel('Stock Price   ','FontSize',12)
zlabel('Inventory Level   ','FontSize',12)
zlim([0,100])
xlim([stock(low),stock(high)])
title('(c)   ')
%view(0,90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)

indent = 20;
X2 = smin + indent*ds:ds:smax-indent*ds;
Y2 = qmin:dq:qmax;

subplot(2,3,1)
h=surf(X2,Y2,marketgrid(:,indent:end-(indent+1),Jt));
set(h,'LineStyle','none')
view(0,90)
%xlim([smin,Js+1])
ylim([0,N])
xlim([smin+indent*ds,smax-indent*ds])
xlabel('Stock Price  ','FontSize',12)
ylabel('Inventory Level  ','FontSize',12)
title('TTM = 0   ','FontSize',10)

subplot(2,3,2)
h=surf(X2,Y2,marketgrid(:,indent:end-(indent+1),Jt-round(Jt/5,0)));
set(h,'LineStyle','none')
view(0,90)
%xlim([smin,Js+1])
ylim([0,N])
xlim([smin+indent*ds,smax-indent*ds])
xlabel('Stock Price  ','FontSize',12)
ylabel('Inventory Level  ','FontSize',12)
title('TTM = 0.2   ','FontSize',10)

subplot(2,3,3)
h4=surf(X2,Y2,marketgrid(:,indent:end-(indent+1),Jt-round(2*Jt/5,0)));
set(h4,'LineStyle','none')
view(0,90)
%xlim([smin,smax])
ylim([0,N])
xlim([smin+indent*ds,smax-indent*ds])
xlabel('Stock Price  ','FontSize',12)
ylabel('Inventory Level  ','FontSize',12)
title('TTM = 0.4   ','FontSize',10)

subplot(2,3,4)
h4=surf(X2,Y2,marketgrid(:,indent:end-(indent+1),Jt-round(3*Jt/5,0)));
set(h4,'LineStyle','none')
view(0,90)
%xlim([smin,smax])
ylim([0,N])
xlim([smin+indent*ds,smax-indent*ds])
xlabel('Stock Price  ','FontSize',12)
ylabel('Inventory Level  ','FontSize',12)
title('TTM = 0.6   ','FontSize',10)

subplot(2,3,5)
h4=surf(X2,Y2,marketgrid(:,indent:end-(indent+1),Jt-round(4*Jt/5,0)));
set(h4,'LineStyle','none')
view(0,90)
%xlim([smin,smax])
ylim([0,N])
xlim([smin+indent*ds,smax-indent*ds])
xlabel('Stock Price  ','FontSize',12)
ylabel('Inventory Level  ','FontSize',12)
title('TTM = 0.8   ','FontSize',10)

subplot(2,3,6)
h4=surf(X2,Y2,marketgrid(:,indent:end-(indent+1),1));
set(h4,'LineStyle','none')
view(0,90)
%xlim([smin,smax])
ylim([0,N])
xlim([smin+indent*ds,smax-indent*ds])
xlabel('Stock Price  ','FontSize',12)
ylabel('Inventory Level  ','FontSize',12)
title('TTM = 1   ','FontSize',10)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
X3 = smin:ds:smax;
Y3 = 0:dq:N;
h5=surf(X3,Y3,grid(21:end-20,:,Jt+1));
set(h5, 'edgecolor','none')
colormap(summer)
freezeColors
%hold on
%h6=surf(X,Y,grid(21:end-20,:,1));
%set(h6, 'edgecolor','none')
%colormap(winter)
freezeColors
%set(h5,'FaceColor',[1 1 1]);
%xlim([smin,smax])
ylim([0,N])
xlim([0,smax])
xlabel('Stock Price  ','FontSize',16)
ylabel('Inventory Level  ','FontSize',16)
zlabel('Terminal Condition (\theta_T)  ','FontSize',16)



















