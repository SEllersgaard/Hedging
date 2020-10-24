%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This programme performs a MC estimate of how well the strategy is doing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global K;
global vol;
global Y;
global kappa;
global lambda;
global g;
global Jt;
global dt;
global ds;
global smin;
global deltaPLUS;           %zeros(Jq+1,Js+1,Jt);   Jt = 1 (means t=0)
global deltaMINUS;          %zeros(Jq+1,Js+1,Jt);   Jt = 1 (means t=0)
global Upper;               %zeros(Jt,Js+1);        Jt = 1 (means t=0)
global Lower;               %zeros(Jt,Js+1);        Jt = 1 (means t=0)
global Mid;                 %zeros(Jt,Js+1);        Jt = 1 (means t=0)
global probability;
global cut;    
global qmin;


initialwealth = 10000; % The initial bank endowment
paths = 10000; %Number of MC simulated paths

Stime = zeros(paths,Jt);                %Stock path

Qtime = zeros(paths,Jt);                %Inventory path
Bank  = zeros(paths,Jt);                %Money account

Qtime2 = zeros(paths,Jt);               %Inventory path
Bank2  = zeros(paths,Jt);               %Money account

Stime(:,1) = K;
s0 = (Stime(1)-smin)/ds +1;
Qtime(:,1) = round(Mid(1,s0));
Bank(:,1)  = initialwealth;

Qtime2(:,1) = round(Mid(1,s0));
Bank2(:,1)  = initialwealth;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a stock path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:1:paths

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
                   Stime(j,t) = Stime(j,t-1) + i*vol;
                   jumpsize(t,1) = i;
                   break
               end
            end
        else
            Stime(j,t) = Stime(j,t-1);
        end

         if tline2(t,1) == 1

            u = rand;
            for i = 1:1:cut
               if u >= probsum(i,1) && u < probsum(i+1,1)
                   Stime(j,t) = Stime(j,t) - i*vol;
                   jumpsize2(t,1) = i;
                   break
               end
            end
         end


    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the evolution of the bank and inventory with limit orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for t = 2:1:Jt
    for j = 1:1:paths
   
        s = (Stime(j,t)-smin)/ds + 1;
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
    

        if (Qtime(j,t-1) -Mid(t,s)) >=1 && Qtime(j,t-1) < Upper(t,s)

            u = rand;
            if u <= lambda*exp(-kappa*dm)*dt
                Bank(j,t) = Bank(j,t-1) + (Stime(j,t)+dm);
                Qtime(j,t) = Qtime(j,t-1)-1;


            else
                Bank(j,t) = Bank(j,t-1);
                Qtime(j,t) = Qtime(j,t-1);
            end

        elseif (Qtime(j,t-1) -Mid(t,s)) >=1 && Qtime(j,t-1) >= Upper(t,s)

            dist = ceil(Qtime(j,t-1)-Upper(t,s));
            Bank(j,t) = Bank(j,t-1) + (Stime(j,t)-Y)*dist;
            Qtime(j,t) = Qtime(j,t-1) - dist;


        elseif (Mid(t,s)-Qtime(j,t-1)) >=1 && Qtime(j,t-1) > Lower(t,s)

            u = rand;
            if u <= lambda*exp(-kappa*dp)*dt
                Bank(j,t) = Bank(j,t-1) - (Stime(j,t)-dp);
                Qtime(j,t) = Qtime(j,t-1)+1;


            else
                Bank(j,t) = Bank(j,t-1);
                Qtime(j,t) = Qtime(j,t-1);
            end

        elseif (Mid(t,s)-Qtime(j,t-1)) >=1 && Qtime(j,t-1) <= Lower(t,s)

            dist = ceil(Lower(t,s)-Qtime(j,t-1));
            Bank(j,t) = Bank(j,t-1) - (Stime(j,t)+Y)*dist;
            Qtime(j,t) = Qtime(j,t-1) + dist;


        else

            Bank(j,t) = Bank(j,t-1);
            Qtime(j,t) = Qtime(j,t-1);

        end
    
    end %end paths
end %end time

% for i = 50:50:paths
%    
%    plot(Stime(i,:));
%    hold on
%     
% end
% hold off
% for i = 50:50:paths
%    
%    plot(Bank(i,:));
%    hold on
%     
% end
% hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the evolution of the bank and inventory with pure market orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for t = 2:1:Jt
   for j=1:1:paths
    s = (Stime(j,t)-smin)/ds + 1;
    
   if Qtime2(j,t-1) - Mid(t,s) >= 1 
        
        dist = ceil(Qtime2(j,t-1)-Mid(t,s));
        Bank2(j,t) = Bank2(j,t-1) + (Stime(j,t)-Y)*dist;
        Qtime2(j,t) = Qtime2(j,t-1) - dist;
       
                
   elseif Mid(t,s) - Qtime2(j,t-1) >= 1  
        
        dist = ceil(Mid(t,s)-Qtime2(j,t-1));
        Bank2(j,t) = Bank2(j,t-1) - (Stime(j,t)+Y)*dist;
        Qtime2(j,t) = Qtime2(j,t-1) + dist;
    
        
    else
        
        Bank2(j,t) = Bank2(j,t-1);
        Qtime2(j,t) = Qtime2(j,t-1);
        
    end
   end  
end



figure(1)
h0 = histogram(Stime(:,end));
h0.Normalization = 'probability';
h0.BinWidth = 2;
%hold on
%x=[K,K];
%y=[0,0.2];
%plot(x,y);
xlabel('Terminal Stock Price  ','FontSize',16)
ylabel('Normalised Frequency  ','FontSize',16)



figure(2)
h1=histogram(Bank(:,end)./initialwealth);
hold on
h2=histogram(Bank2(:,end)./initialwealth);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
xlabel('B_T/B_0  ','FontSize',16)
ylabel('Normalised Frequency  ','FontSize',16)
h_legend=legend('With Limit Orders   ','Without Limit Orders   ');
set(h_legend,'FontSize',14);
% figure(3)
% for j = 1:1:paths
% plot(Bank(j,:))
% hold on
% end
% hold off
% 
% figure(4)
% for j = 1:1:paths
% plot(Bank2(j,:))
% hold on
% end
% hold off

% figure(5)
% for j = 1:1:paths
% plot(Qtime(j,:))
% hold on
% end
% hold off
% figure(6)
% for j = 1:1:paths
% plot(Qtime2(j,:))
% hold on
% end
% hold off

C = cell(8,3);
C{1,1} = '';
C{1,2} = 'Limit';
C{1,3} = 'Market';
C{2,1} = 'Mean';
C{2,2} = mean(Bank(:,end)./initialwealth);
C{2,3} = mean(Bank2(:,end)./initialwealth);
C{3,1} = 'Std';
C{3,2} = std(Bank(:,end)./initialwealth);
C{3,3} = std(Bank2(:,end)./initialwealth);
C{4,1} = '1% Quantile';
C{4,2} = quantile(Bank(:,end)./initialwealth,0.01);
C{4,3} = quantile(Bank2(:,end)./initialwealth,0.01);
C{5,1} = '25% Quantile';
C{5,2} = quantile(Bank(:,end)./initialwealth,0.25);
C{5,3} = quantile(Bank2(:,end)./initialwealth,0.25);
C{6,1} = '50% Quantile';
C{6,2} = quantile(Bank(:,end)./initialwealth,0.5);
C{6,3} = quantile(Bank2(:,end)./initialwealth,0.5);
C{7,1} = '75% Quantile';
C{7,2} = quantile(Bank(:,end)./initialwealth,0.75);
C{7,3} = quantile(Bank2(:,end)./initialwealth,0.75);
C{8,1} = '99% Quantile';
C{8,2} = quantile(Bank(:,end)./initialwealth,0.99);
C{8,3} = quantile(Bank2(:,end)./initialwealth,0.99);
C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include the terminal condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wealth  = zeros(paths,1);
Wealth2  = zeros(paths,1);



for j = 1:1:paths
   
    s = (Stime(j,end)-smin)/ds + 1;
    S = Stime(j,end);
    B = Bank(j,end);
    Q = Qtime(j,end);
    
    if Q >= Mid(end,s) && S >= K
        Wealth(j,1) = B + N*K + (Q-N)*(S-Y);
    elseif  Q >= Mid(end,s) && S < K
        Wealth(j,1) = B + Q*(S-Y);
    elseif  Q < Mid(end,s) && S >= K
        Wealth(j,1) = B + N*K - (N-Q)*(S+Y);
    elseif  Q < Mid(end,s) && S < K
        Wealth(j,1) = B + Q*(S-Y);
    end
    
    B2 = Bank2(j,end);
    Q2 = Qtime2(j,end);
    
    if Q2 >= Mid(end,s) && S >= K
        Wealth2(j,1) = B2 + N*K + (Q2-N)*(S-Y);
    elseif  Q2 >= Mid(end,s) && S < K
        Wealth2(j,1) = B2 + Q2*(S-Y);
    elseif  Q2 < Mid(end,s) && S >= K
        Wealth2(j,1) = B2 + N*K - (N-Q2)*(S+Y);
    elseif  Q2 < Mid(end,s) && S < K
        Wealth2(j,1) = B2 + Q2*(S-Y);
    end

end

callprice = sqrt(2*g)*(sqrt(2-probability)/probability)*vol*normpdf(0,0,1);
wstock = round(Mid(1,s0))*K - 100*callprice;

figure(3)
h3=histogram(Wealth(:,end)./(initialwealth+wstock));
hold on
h4=histogram(Wealth2(:,end)./(initialwealth+wstock));
h3.Normalization = 'probability';
h3.BinWidth = 0.005;
h4.Normalization = 'probability';
h4.BinWidth = 0.005;
xlabel('W_T/B_0  ','FontSize',16)
ylabel('Normalised Frequency  ','FontSize',16)
h_legend=legend('With Limit Orders   ','Without Limit Orders   ');
set(h_legend,'FontSize',14);

D = cell(8,3);
D{1,1} = '';
D{1,2} = 'Limit';
D{1,3} = 'Market';
D{2,1} = 'Mean';
D{2,2} = mean(Wealth(:,end)./(initialwealth+wstock));
D{2,3} = mean(Wealth2(:,end)./(initialwealth+wstock));
D{3,1} = 'Std';
D{3,2} = std(Wealth(:,end)./(initialwealth+wstock));
D{3,3} = std(Wealth2(:,end)./(initialwealth+wstock));
D{4,1} = '1% Quantile';
D{4,2} = quantile(Wealth(:,end)./(initialwealth+wstock),0.01);
D{4,3} = quantile(Wealth2(:,end)./(initialwealth+wstock),0.01);
D{5,1} = '25% Quantile';
D{5,2} = quantile(Wealth(:,end)./(initialwealth+wstock),0.25);
D{5,3} = quantile(Wealth2(:,end)./(initialwealth+wstock),0.25);
D{6,1} = '50% Quantile';
D{6,2} = quantile(Wealth(:,end)./(initialwealth+wstock),0.5);
D{6,3} = quantile(Wealth2(:,end)./(initialwealth+wstock),0.5);
D{7,1} = '75% Quantile';
D{7,2} = quantile(Wealth(:,end)./(initialwealth+wstock),0.75);
D{7,3} = quantile(Wealth2(:,end)./(initialwealth+wstock),0.75);
D{8,1} = '99% Quantile';
D{8,2} = quantile(Wealth(:,end)./(initialwealth+wstock),0.99);
D{8,3} = quantile(Wealth2(:,end)./(initialwealth+wstock),0.99);
D




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hFig = figure(4);

subplot(1,3,1)
%h0 = histogram(Stime(:,end),'Normalization','probability');
%h0.Normalization = 'probability';
%h0.BinWidth = 2;
h = histfit(Stime(:,end),40);
h(1).FaceColor = [0.5 0.5 0.5];
h(2).Color = [0 0 0];
xlabel('S_T  ','FontSize',12)
ylabel('Frequency  ','FontSize',12)
title('(a)   ','FontSize',10)
% mu = mean(Stime(:,end));
% sigma = sqrt(2*g)*(sqrt(2-probability)/probability)*vol;  %std(Stime(:,end))
% hold on
% y = 0:0.1:100;
% f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% plot(y,f,'LineWidth',1.5)

subplot(1,3,2)
h1=histogram(Bank(:,end)./initialwealth);
hold on
h2=histogram(Bank2(:,end)./initialwealth);
%h1.Normalization = 'probability';
h1.BinWidth = 0.02;
h1.FaceColor = [0.5 0.5 0.5];
%h2.Normalization = 'probability';
h2.BinWidth = 0.02;
h2.FaceColor = [0.1 0.1 0.1];
xlabel('B_T/B_0  ','FontSize',12)
ylabel('Frequency  ','FontSize',12)
h_legend=legend('With Limit Orders   ','Without Limit Orders   ');
%set(h_legend,'FontSize',14);
title('(b)   ','FontSize',10)

subplot(1,3,3)

h3=histogram(Wealth(:,end)./(initialwealth+wstock));
hold on
h4=histogram(Wealth2(:,end)./(initialwealth+wstock));
%h3.Normalization = 'probability';
h3.BinWidth = 0.005;
h3.FaceColor = [0.5 0.5 0.5];
%h4.Normalization = 'probability';
h4.BinWidth = 0.005;
h4.FaceColor = [0.1 0.1 0.1];
xlabel('\Pi_T/\Pi_0^B  ','FontSize',12)
ylabel('Frequency  ','FontSize',12)
h_legend=legend('With Limit Orders   ','Without Limit Orders   ');
title('(c)   ','FontSize',10)


set(hFig, 'Position', [50 50 1200 300])



