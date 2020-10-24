%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This programme simulates a stock path and computes the hedge strategy
% based on the optimal boudaries. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% global K;
% global vol;
% global T;
% global Y;
% global kappa;
% global lambda;
% global g;
% global Jt;
% global dt;
% global ds;
% global smin;
% global smax;
% global deltaPLUS;           %zeros(Jq+1,Js+1,Jt);   Jt = 1 (means t=0)
% global deltaMINUS;          %zeros(Jq+1,Js+1,Jt);   Jt = 1 (means t=0)
% global Upper;               %zeros(Jt,Js+1);        Jt = 1 (means t=0)
% global Lower;               %zeros(Jt,Js+1);        Jt = 1 (means t=0)
% global Mid;                 %zeros(Jt,Js+1);        Jt = 1 (means t=0)
% 
% 
% Stime = zeros(Jt,1);                %Stock path
% Qtime = zeros(Jt,1);                %Inventory path
% Bank  = zeros(Jt,1);                %Money account
% limit = zeros(Jt,1);                %When are limit orders successful 
% market = zeros(Jt,1);               %When are market orders being placed
% market2 = zeros(Jt,1);
% 
% Qtime2 = zeros(Jt,1);                %Inventory path
% Bank2  = zeros(Jt,1);                %Money account
% 
% Stime(1) = K;
% s0 = (Stime(1)-smin)/ds +1;
% Qtime(1) = round(Mid(1,s0));
% Bank(1)  = 10000;
% Qtime2(1) = round(Mid(1,s0));
% Bank2(1)  = 10000;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Generate a stock path
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% for t = 2:1:Jt
%    
%     u1 = rand;
%     u2 = rand;
%     
%     if u1 <= g*dt   
%         Stime(t) = Stime(t-1) + vol;
%     else
%         Stime(t) = Stime(t-1);
%     end    
%     if u2 <= g*dt
%         Stime(t) = Stime(t) - vol;
%     end
% 
% end

global K;
global vol;
global T;
global Y;
global kappa;
global lambda;
global g;
global Jt;
global dt;
global ds;
global smin;
global smax;
global deltaPLUS;           %zeros(Jq+1,Js+1,Jt);   Jt = 1 (means t=0)
global deltaMINUS;          %zeros(Jq+1,Js+1,Jt);   Jt = 1 (means t=0)
global Upper;               %zeros(Jt,Js+1);        Jt = 1 (means t=0)
global Lower;               %zeros(Jt,Js+1);        Jt = 1 (means t=0)
global Mid;                 %zeros(Jt,Js+1);        Jt = 1 (means t=0)
global probability;
global cut;
global qmin;

Stime = zeros(Jt,1);                %Stock path
Qtime = zeros(Jt,1);                %Inventory path
Bank  = zeros(Jt,1);                %Money account
limit = zeros(Jt,1);                %When are limit orders successful 
market = zeros(Jt,1);               %When are market orders being placed
market2 = zeros(Jt,1);

Qtime2 = zeros(Jt,1);                %Inventory path
Bank2  = zeros(Jt,1);                %Money account

Stime(1) = K;
s0 = (Stime(1)-smin)/ds +1;
Qtime(1) = round(Mid(1,s0));
Bank(1)  = 10000;
Qtime2(1) = round(Mid(1,s0));
Bank2(1)  = 10000;

delp  = zeros(Jt,1);                
delm  = zeros(Jt,1);





tsum = 0;
tline = zeros(Jt,1);
tsum2 = 0;
tline2 = zeros(Jt,1);

while 1
   
    u1 = rand;
    decay1 = -(1/g)*log(1-u1);
    tsum = tsum + decay1;
    if tsum >= 1
        break
    end
    stepdecay = max(round(tsum/dt),1);
    tline(stepdecay,1) = 1;
    
end

while 1
   
    u2 = rand;
    decay2 = -(1/g)*log(1-u2);
    tsum2 = tsum2 + decay2;
    if tsum2 >= 1
        break
    end
    stepdecay2 = max(round(tsum2/dt),1);
    tline2(stepdecay2,1) = 1;
 
end

probsum = zeros(cut+1,1);
prs = 0;
for i = 2:1:cut
    prs = prs + (1-probability).^(i-2)*probability;
    probsum(i,1) = prs;
end
probsum(end,1) = 1;


jumpsize = zeros(Jt,1);
jumpsize2 = zeros(Jt,1);


for t = 2:1:Jt
   
    if tline(t,1) == 1
        
        u = rand;
        for i = 1:1:cut
           if u >= probsum(i,1) && u < probsum(i+1,1)
               Stime(t) = Stime(t-1) + i*vol;
               jumpsize(t,1) = i;
               break
           end
        end
    else
        Stime(t) = Stime(t-1);
    end
    
     if tline2(t,1) == 1
        
        u = rand;
        for i = 1:1:cut
           if u >= probsum(i,1) && u < probsum(i+1,1)
               Stime(t) = Stime(t) - i*vol;
               jumpsize2(t,1) = i;
               break
           end
        end
     end
   
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the evolution of the bank and inventory with limit orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

limitsell = 0;
limitbuy = 0;

for t = 2:1:Jt
   
    s = (Stime(t)-smin)/ds + 1;
    dm = deltaMINUS(abs(qmin)+Qtime(t-1)+1,s,t);
    dp = deltaPLUS(abs(qmin)+Qtime(t-1)+1,s,t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IN CASE A ZERO DELTA IS SELECTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if deltaMINUS(abs(qmin)+Qtime(t-1)+1,s,t) == 0
        dm = 0;
        count1 = 0;
        while dm == 0
            count1 = count1 + 1;
            dm = deltaMINUS(abs(qmin)+Qtime(t-1)+1+count1,s,t);
        end    
    end

    if deltaPLUS(abs(qmin)+Qtime(t-1)+1,s,t) == 0
        dp = 0;
        count2 = 0;
        while dp == 0
            count2 = count2 - 1;
            dp = deltaPLUS(abs(qmin)+Qtime(t-1)+1+count2,s,t);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if (Qtime(t-1) -Mid(t,s)) >=1 && Qtime(t-1) < Upper(t,s)
        
        u = rand;
        if u <= lambda*exp(-kappa*dm)*dt
            Bank(t) = Bank(t-1) + (Stime(t)+dm);
            delm(t) = dm;
            Qtime(t) = Qtime(t-1)-1;
            limitsell = limitsell + 1;
            limit(t) = -1;
            %[-1,abs(qmin)+Qtime(t-1)+1,s,t,deltaMINUS(abs(qmin)+Qtime(t-1)+1,s,t),dm]
            
        else
            Bank(t) = Bank(t-1);
            Qtime(t) = Qtime(t-1);
        end
        
    elseif (Qtime(t-1) -Mid(t,s)) >=1 && Qtime(t-1) >= Upper(t,s)
        
        dist = ceil(Qtime(t-1)-Upper(t,s));
        Bank(t) = Bank(t-1) + (Stime(t)-Y)*dist;
        Qtime(t) = Qtime(t-1) - dist;
        market(t) = -dist;
        
    elseif (Mid(t,s)-Qtime(t-1)) >=1 && Qtime(t-1) > Lower(t,s)
        
        u = rand;
        if u <= lambda*exp(-kappa*dp)*dt
            Bank(t) = Bank(t-1) - (Stime(t)-dp);
            delp(t) = dp;
            Qtime(t) = Qtime(t-1)+1;
            limitbuy = limitbuy + 1;
            limit(t) = 1;
            %[1,abs(qmin)+Qtime(t-1)+1,s,t,deltaPLUS(abs(qmin)+Qtime(t-1)+1,s,t),dp]
        else
            Bank(t) = Bank(t-1);
            Qtime(t) = Qtime(t-1);
        end
        
    elseif (Mid(t,s)-Qtime(t-1)) >=1 && Qtime(t-1) <= Lower(t,s)
        
        dist = ceil(Lower(t,s)-Qtime(t-1));
        Bank(t) = Bank(t-1) - (Stime(t)+Y)*dist;
        Qtime(t) = Qtime(t-1) + dist;
        market(t) = dist;
        
    else
        
        Bank(t) = Bank(t-1);
        Qtime(t) = Qtime(t-1);
        
    end
      
end


time = 0:dt:(T-dt);

% figure(1)
% 
% subplot(2,2,1)
% plot(time,Stime,'r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Stock Price  ','FontSize',16)
% %title('Stock Path Evolution','FontSize',16)
% 
% subplot(2,2,2)
% plot(time,Qtime,'r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Inventory Level  ','FontSize',16)
% %title('Inventory Level Evolution','FontSize',16)
% 
% subplot(2,2,3)
% plot(time,Bank,'r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Bank Account  ','FontSize',16)
% %title('Stock Path Evolution','FontSize',16)
% 
% subplot(2,2,4)
% stem(time,limit,'.k');
% hold on
% stem(time,market,'.r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Order Executions  ','FontSize',16)
% legend('Limit Order  ','Market Order  ')
% 
% p = mtit('Hedging with Limit and Market Orders  ','FontSize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the evolution of the bank and inventory with pure market orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t = 2:1:Jt
   
    s = (Stime(t)-smin)/ds + 1;
    
   if Qtime2(t-1) - Mid(t,s) >= 1 
        
        dist = ceil(Qtime2(t-1)-Mid(t,s));
        Bank2(t) = Bank2(t-1) + (Stime(t)-Y)*dist;
        Qtime2(t) = Qtime2(t-1) - dist;
        market2(t) = -dist;
                
   elseif Mid(t,s) - Qtime2(t-1) >= 1  
        
        dist = ceil(Mid(t,s)-Qtime2(t-1));
        Bank2(t) = Bank2(t-1) - (Stime(t)+Y)*dist;
        Qtime2(t) = Qtime2(t-1) + dist;
        market2(t) = dist;
        
    else
        
        Bank2(t) = Bank2(t-1);
        Qtime2(t) = Qtime2(t-1);
        
    end
      
end

% 
% figure(2)
% 
% subplot(2,2,1)
% plot(time,Stime,'r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Stock Price  ','FontSize',16)
% 
% subplot(2,2,2)
% plot(time,Qtime2,'r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Inventory Level  ','FontSize',16)
% %title('Inventory Level Evolution','FontSize',16)
% 
% subplot(2,2,3)
% plot(time,Bank2,'r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Bank Account  ','FontSize',16)
% 
% subplot(2,2,4)
% stem(time,market2,'.r');
% xlabel('Time [Years]  ','FontSize',16)
% ylabel('Order Executions  ','FontSize',16)
% 
% q = mtit('Hedging with Market Orders Only  ','FontSize',16);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hFig = figure(3);

clp = [0.6 0.6 0.6];

subplot(2,3,1)
plot(time,Stime,'k');
xlabel('Time [Years]  ','FontSize',12)
ylabel('Stock Price (S)  ','FontSize',12)
title('(a)   ','FontSize',10)

subplot(2,3,2)
plot(time,Qtime2,'Color',clp);
hold on
plot(time,Qtime,'k')
xlabel('Time [Years]  ','FontSize',12)
ylabel('Inventory Level (Q)  ','FontSize',12)
title('(b)   ','FontSize',10)
legend('Market Orders  ','Lim + Market Orders  ')

subplot(2,3,3)
plot(time,Bank2,'Color',clp);
hold on
plot(time,Bank,'k')
xlabel('Time [Years]  ','FontSize',12)
ylabel('Bank Account (B)  ','FontSize',12)
title('(c)   ','FontSize',10)
legend('Market Orders  ','Lim + Market Orders  ')

subplot(2,3,4)
stem(time,market2,'.','Color',clp);
xlabel('Time [Years]  ','FontSize',12)
ylabel('Order Executions  ','FontSize',12)
title('(d)   ','FontSize',10)

subplot(2,3,5)
market3 = market;
limit2 = limit;
market3(market3 == 0) = NaN;
limit2(limit2 == 0) = NaN;
stem(time,market3,'.','Color',clp,'ShowBaseLine','off');
hold on
stem(time,limit2,'xk','ShowBaseLine','off');
hold on
rl = refline(0,0);
rl.Color= 'k';
xlabel('Time [Years]  ','FontSize',12)
ylabel('Order Executions  ','FontSize',12)
title('(e)   ','FontSize',10)
legend('Market Order  ','Limit Order  ')

subplot(2,3,6)
delm2 = delm;
delp2 = delp;
delm2(delm2==0) = NaN;
delp2(delp2==0) = NaN;
stem(time,delm2,'xk');
hold on
stem(time,delp2,'.','Color',clp);
xlabel('Time [Years]  ','FontSize',12)
ylabel('Spread (\delta)  ','FontSize',12)
title('(f)   ','FontSize',10)
legend('\delta^-  ','\delta^+  ')

set(hFig, 'Position', [50 50 900 500])

% subplot(3,3,[7,8,9])
% plot(time,Stime,'k');
% xlabel('Time [Years]  ','FontSize',12)
% ylabel('Stock Price (S)  ','FontSize',12)
% title('(a)   ','FontSize',10)
% hold on
% plot(time,Stime-2.5,'k');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make a 3D stem plot of the inventory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InventoryGrid = zeros(Jt,Js+1);
InventoryGrid2 = zeros(Jt,Js+1);

for t = 1:1:Jt
   
    coord = round((Stime(t)+ds)/ds);
    
    InventoryGrid(t,coord) = Qtime(t);
    
    InventoryGrid2(t,coord) = Qtime2(t);
    
    %%%%FOR VISUAL REASONS ONLY (SEE BELOW)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if InventoryGrid(t,coord) == 0
        InventoryGrid(t,coord) = 1;
    end
    if InventoryGrid2(t,coord) == 0
        InventoryGrid2(t,coord) = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
InventoryGrid(InventoryGrid == 0) = NaN;   % Do this - otherwise stem3 will display zero points
InventoryGrid2(InventoryGrid2 == 0) = NaN;

 
 
 time2 = T-dt:-dt:0;
 stock2 = smin:ds:smax;
 
% figure(3)
% 
% h = surf(stock2,time2,Mid);
% set(h,'LineStyle','none')
% set(h,'FaceColor',[1 1 0],'FaceAlpha',0.5)
% hold on
% h = stem3(stock2,time2,InventoryGrid);
% 
% ylabel('TTM   ','FontSize',16)
% xlabel('Stock price   ','FontSize',16)
% zlabel('Market Order Boundary [# of Stocks]  ','FontSize',16)
% zlim([0,100])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the same but with fewer data points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



indent = 20;
X2 = smin + indent*ds:ds:smax-indent*ds;
Y2 = qmin:dq:qmax;
stock2 = X2;



IG = zeros(Jt/10+1,length(X2));
IG2 = zeros(Jt/10+1,length(X2));
IM = zeros(Jt/10+1,length(X2));
IU = zeros(Jt/10+1,length(X2));
IL = zeros(Jt/10+1,length(X2));

count = 0;
IG(1,:) = InventoryGrid(1,indent:end-(indent+1));
IG2(1,:) = InventoryGrid2(1,indent:end-(indent+1));
IM(1,:) = Mid(1,indent:end-(indent+1));
for t = 10:10:Jt
   count = count +1; 
   IG(count,:) = InventoryGrid(t,indent:end-(indent+1));
   IG2(count,:) = InventoryGrid2(t,indent:end-(indent+1));
   IM(count,:) = Mid(t,indent:end-(indent+1));
   IU(count,:) = Upper(t,indent:end-(indent+1));
   IL(count,:) = Lower(t,indent:end-(indent+1));
end
IM(end,:) = IM(end-1,:);
IU(end,:) = IU(end-1,:);
IL(end,:) = IL(end-1,:);
IG(IG == 0) = NaN;
IG2(IG2 == 0) = NaN;


time3 = [T-dt:-10*dt:0,0];

%%

ucl = [0.7 0.7 0.7];
lcl = [0.4 0.4 0.4];
mcl = [1 1 1];

lop = 0.7;
uop = 0.3;

dxx = 3;
dyy = 3;

hFig = figure(4);
subplot(1,3,1)
h = surf(stock2(1:dxx:end),time3(1:dyy:end),IM(1:dyy:end,1:dxx:end));
%set(h,'LineStyle','none')
set(h,'FaceAlpha',0)
set(h,'FaceColor',mcl)
hold on
h2 = stem3(stock2,time3,IG,'filled');
h2.Color = 'r';
h2.MarkerFaceColor = 'k';
ylabel('TTM   ','FontSize',12)
xlabel('Stock (S)   ','FontSize',12)
zlabel('Inventory Level (Q)  ','FontSize',12)
zlim([0,100])
xlim([smin+indent*ds,smax-indent*ds])
title('(a)  ','FontSize',12)
freezeColors 

hold on 
h99 = surf(stock2,time3,IU);
set(h99,'LineStyle','none')
set(h99,'FaceAlpha',uop)
%colormap(bone)
%set(h99,'FaceColor',[53/255 42/255 134/255])
set(h99,'FaceColor',ucl)
freezeColors 

hold on 
h999 = surf(stock2,time3,IL);
set(h999,'LineStyle','none')
set(h999,'FaceAlpha',lop)
%set(h999,'FaceColor',[166/255 191/255 106/255])
set(h999,'FaceColor',lcl)
%colormap(flipud(colormap(summer)))
freezeColors 


subplot(1,3,2)
h = surf(stock2(1:dxx:end),time3(1:dyy:end),IM(1:dyy:end,1:dxx:end));
%set(h,'LineStyle','none')
set(h,'FaceColor',mcl)
set(h,'FaceAlpha',0)
%colormap(autumn) 
hold on
h2 = stem3(stock2,time3,IG,'filled');
h2.Color = 'r';
h2.MarkerFaceColor = 'k';
ylabel('TTM   ','FontSize',12)
xlabel('Stock Price (S)   ','FontSize',12)
zlabel('Inventory Level (Q)  ','FontSize',12)
zlim([0,100])
xlim([smin+indent*ds,smax-indent*ds])
title('(b)  ','FontSize',12)
freezeColors 

hold on 
h99 = surf(stock2,time3,IU);
set(h99,'LineStyle','none')
set(h99,'FaceAlpha',uop)
%colormap(bone)
%set(h99,'FaceColor',[53/255 42/255 134/255])
set(h99,'FaceColor',ucl)
freezeColors 

hold on 
h999 = surf(stock2,time3,IL);
set(h999,'LineStyle','none')
set(h999,'FaceAlpha',lop)
%set(h999,'FaceColor',[166/255 191/255 106/255])
set(h999,'FaceColor',lcl)
%colormap(flipud(colormap(summer)))
freezeColors 
view(0,90)


%%%%%%

subplot(1,3,3)
h = surf(stock2(1:dxx:end),time3(1:dyy:end),IM(1:dyy:end,1:dxx:end));
%set(h,'LineStyle','none')
set(h,'FaceColor',mcl)
set(h,'FaceAlpha',0)
%colormap(autumn) 
hold on
h2 = stem3(stock2,time3,IG,'filled');
h2.Color = 'r';
h2.MarkerFaceColor = 'k';
ylabel('TTM   ','FontSize',12)
xlabel('Stock Price (S)   ','FontSize',12)
zlabel('Inventory Level (Q)  ','FontSize',12)
zlim([0,100])
xlim([smin+indent*ds,smax-indent*ds])
title('(c)  ','FontSize',12)
freezeColors 

hold on 
h99 = surf(stock2,time3,IU);
set(h99,'LineStyle','none')
set(h99,'FaceAlpha',uop)
%colormap(bone)
%set(h99,'FaceColor',[53/255 42/255 134/255])
set(h99,'FaceColor',ucl)
freezeColors 

hold on 
h999 = surf(stock2,time3,IL);
set(h999,'LineStyle','none')
set(h999,'FaceAlpha',lop)
%set(h999,'FaceColor',[166/255 191/255 106/255])
set(h999,'FaceColor',lcl)
%colormap(flipud(colormap(summer)))
freezeColors 
view(180,0)

set(hFig, 'Position', [50 50 1200 300])

%%%%PLOT SURFACE CONTOURS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% plot3(stock2,ones(1,length(IU)),IU(1,:),'k','LineWidth',1)
% hold on
% plot3(stock2,ones(1,length(IL)),IL(1,:),':k','LineWidth',1)
% hold on
% plot3(stock2,0.75*T*ones(1,length(IU)),IU(25,:),'k','LineWidth',1)
% hold on
% plot3(stock2,0.75*T*ones(1,length(IL)),IL(25,:),':k','LineWidth',1)
% hold on
% plot3(stock2,0.5*T*ones(1,length(IU)),IU(50,:),'k','LineWidth',1)
% hold on
% plot3(stock2,0.5*T*ones(1,length(IL)),IL(50,:),':k','LineWidth',1)
% hold on
% plot3(stock2,0.25*T*ones(1,length(IU)),IU(75,:),'k','LineWidth',1)
% hold on
% plot3(stock2,0.25*T*ones(1,length(IL)),IL(75,:),':k','LineWidth',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(6)
% h = surf(stock2,time3,IM);
% set(h,'LineStyle','none')
% set(h,'FaceAlpha',0.5)
% colormap(autumn) 
% hold on
% h2 = stem3(stock2,time3,IG2,'filled');
% h2.Color = 'r';
% h2.MarkerFaceColor = 'k';
% ylabel('TTM   ','FontSize',16)
% xlabel('Stock price   ','FontSize',16)
% zlabel('Market Order Boundary [# of Stocks]  ','FontSize',16)
% zlim([0,100])
% xlim([smin+indent*ds,smax-indent*ds])
% title('Hedging with Market Orders Only  ','FontSize',12)










% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Do the same but with fewer data points
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% IG = zeros(Jt/10+1,Js+1);
% IG2 = zeros(Jt/10+1,Js+1);
% IM = zeros(Jt/10+1,Js+1);
% IU = zeros(Jt/10+1,Js+1);
% IL = zeros(Jt/10+1,Js+1);
% 
% count = 0;
% IG(1,:) = InventoryGrid(1,:);
% IG2(1,:) = InventoryGrid2(1,:);
% IM(1,:) = Mid(1,:);
% for t = 10:10:Jt
%    count = count +1; 
%    IG(count,:) = InventoryGrid(t,:);
%    IG2(count,:) = InventoryGrid2(t,:);
%    IM(count,:) = Mid(t,:);
%    IU(count,:) = Upper(t,:);
%    IL(count,:) = Lower(t,:);
% end
% IM(end,:) = IM(end-1,:);
% IU(end,:) = IU(end-1,:);
% IL(end,:) = IL(end-1,:);
% IG(IG == 0) = NaN;
% IG2(IG2 == 0) = NaN;
% 
% 
% time3 = [T-dt:-10*dt:0,0];
% 
% figure(3)
% h = surf(stock2,time3,IM);
% set(h,'LineStyle','none')
% set(h,'FaceAlpha',0.5)
% colormap(autumn) 
% hold on
% h2 = stem3(stock2,time3,IG,'filled');
% h2.Color = 'r';
% h2.MarkerFaceColor = 'k';
% ylabel('TTM   ','FontSize',16)
% xlabel('Stock price   ','FontSize',16)
% zlabel('Market Order Boundary [# of Stocks]  ','FontSize',16)
% zlim([0,100])
% title('Hedging with Limit and Market Orders  ','FontSize',16)
% %%%%PLOT SURFACE CONTOURS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% plot3(stock2,ones(1,length(IU)),IU(1,:),'k','LineWidth',1)
% hold on
% plot3(stock2,ones(1,length(IL)),IL(1,:),':k','LineWidth',1)
% hold on
% plot3(stock2,0.75*T*ones(1,length(IU)),IU(25,:),'k','LineWidth',1)
% hold on
% plot3(stock2,0.75*T*ones(1,length(IL)),IL(25,:),':k','LineWidth',1)
% hold on
% plot3(stock2,0.5*T*ones(1,length(IU)),IU(50,:),'k','LineWidth',1)
% hold on
% plot3(stock2,0.5*T*ones(1,length(IL)),IL(50,:),':k','LineWidth',1)
% hold on
% plot3(stock2,0.25*T*ones(1,length(IU)),IU(75,:),'k','LineWidth',1)
% hold on
% plot3(stock2,0.25*T*ones(1,length(IL)),IL(75,:),':k','LineWidth',1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% figure(4)
% h = surf(stock2,time3,IM);
% set(h,'LineStyle','none')
% set(h,'FaceAlpha',0.5)
% colormap(autumn) 
% hold on
% h2 = stem3(stock2,time3,IG2,'filled');
% h2.Color = 'r';
% h2.MarkerFaceColor = 'k';
% ylabel('TTM   ','FontSize',16)
% xlabel('Stock price   ','FontSize',16)
% zlabel('Market Order Boundary [# of Stocks]  ','FontSize',16)
% zlim([0,100])
% title('Hedging with Market Orders Only  ','FontSize',16)
% 
% 
% % figure(5)
% % h = surf(stock2,time3,IM);
% % set(h,'LineStyle','none')
% % set(h,'FaceAlpha',0.4)
% % colormap(autumn) 
% % freezeColors
% % hold on
% % h = surf(stock2,time3,IU);
% % set(h,'LineStyle','none')
% % set(h,'FaceAlpha',0.3)
% % colormap(winter)
% % freezeColors
% % hold on
% % h = surf(stock2,time3,IL);
% % set(h,'LineStyle','none')
% % set(h,'FaceAlpha',0.5)
% % colormap(flipud(colormap(summer)))
% % hold on
% % 
% % h2 = stem3(stock2,time3,IG,'filled');
% % h2.Color = 'r';
% % h2.MarkerFaceColor = 'k';
% % ylabel('TTM   ','FontSize',16)
% % xlabel('Stock price   ','FontSize',16)
% % zlabel('Market Order Boundary [# of Stocks]  ','FontSize',16)
% % zlim([0,100])
% % title('Hedging with Limit and Market Orders  ','FontSize',16)
% % 
% % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stat


countmarket = 0;
countmarket2 = 0;
countlimit = 0;
countdelm = 0;
countdelp = 0;

for j=1:1:Jt
   
    if market2(j,1) ~= 0
        countmarket2 = countmarket2 +1;
    end
    if market(j,1) ~= 0
        countmarket = countmarket +1;
    end
    if limit(j,1) ~= 0
        countlimit = countlimit +1;
    end
    if delm(j,1) ~= 0
        countdelm = countdelm +1;
    end
    if delp(j,1) ~= 0
        countdelp = countdelp +1;
    end
    
end

[countmarket2,countmarket,countlimit,(sum(delm)+sum(delp))/(countdelm+countdelp)]


DeltaQ = Qtime(2:end,1)-Qtime(1:end-1,1);
DeltaQ2 = Qtime2(2:end,1)-Qtime2(1:end-1,1);
DeltaBank = Bank(2:end,1)-Bank(1:end-1,1);
DeltaBank2 = Bank2(2:end,1)-Bank2(1:end-1,1);

QVQ = mean(DeltaQ.^2);
QVQ2 = mean(DeltaQ2.^2);
QVB = mean(DeltaBank.^2);
QVB2 = mean(DeltaBank2.^2);

[QVQ,QVQ2]

[QVB,QVB2]








