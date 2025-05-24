% 0515-Make all var LN and continuous
% 0508-denotes are from E:\OneDrive - University at Buffalo\2-Buffalo\1-Course\1-Kang Research\5-Send to Kang\wk13-0421-meeting\Cont_Extended_Abstract_0424
% 0424-operational formulation explain look at E:\OneDrive - University at Buffalo\2-Buffalo\1-Course\1-Kang Research\3.Zhiheng Data\2020.04.05 Formulation and Code "Formulation.ppt"
% 0508-hide c17 that limit Uijt in operational, otherwise no solution.[this action is correct, this c should put in upper for next stage]

%0-read all test case 12 [done by wk2-02/07/20fri]

%0-start recording running time
tic

%0.1-read and store test case 12
clear;clc;close all;

% 111820-Iteration Count
Iter=1;

%-----for small range before whole dataset, set i,j,t upper bound------
i_up=3; %should be 20
j_up=3; %should be 20
t_up=4; %should be 40
%--------

d= zeros(20,20,40);%create 20*20*40 matrix
% 101520-read different cijt
c= zeros(20,20,40);
% 101520-pricing range for every i,j,t
P= zeros(20,20,40);
% 111820-Avg P from t=1 to 40
% AvgP = zeros(40);
% 111820-Avg cijt from t=1 to 40
% AvgC = zeros(40);

% read dijt
for t=1:40 %40 time periods
    fileName=['Period',num2str(t-1),'.txt']; %txt: Period0 to 39
    read=load(fileName); %load that txt
    d(:,:,t)=read; %assign txt to dijt
end

% 101520-read cijt and keep standard $19.62/hrs
for t1=1:40 %40 time periods
    fileName1=['cijtPeriod',num2str(t1-1),'.txt']; 
    read1=load(fileName1); %load that txt
    c(:,:,t1)=read1./60.*19.62; %convert traveling time (min/trip) to ($/trip)
    travt(:,:,t1)=read1./60; %travel time: hrs
end

% 061421-std.P=24$/hr, initial pricing exactly to be 24$/hr (33% above std. cost=18$/hr) - all be same

for i2=1:20
  for j2=1:20
    for t2=1:40
      P_inipertrip(i2,j2,t2) = 24.*travt(i2,j2,t2); %all std. pricing = $24/hrs (this P = $/hrs)
%       %{
      if P_inipertrip(i2,j2,t2) < 4
        P_inipertrip(i2,j2,t2) = 2;
      elseif P_inipertrip(i2,j2,t2)>=4 && P_inipertrip(i2,j2,t2)<8
        P_inipertrip(i2,j2,t2)=6;
      elseif P_inipertrip(i2,j2,t2)>=8 && P_inipertrip(i2,j2,t2)<12
        P_inipertrip(i2,j2,t2)=10;  
      elseif P_inipertrip(i2,j2,t2)>=12 && P_inipertrip(i2,j2,t2)<16
        P_inipertrip(i2,j2,t2)=14;
      elseif P_inipertrip(i2,j2,t2)>=16 && P_inipertrip(i2,j2,t2)<20
        P_inipertrip(i2,j2,t2)=18;
      elseif P_inipertrip(i2,j2,t2)>=20 && P_inipertrip(i2,j2,t2)<24
        P_inipertrip(i2,j2,t2)=22;
      elseif P_inipertrip(i2,j2,t2)>=24 && P_inipertrip(i2,j2,t2)<28
        P_inipertrip(i2,j2,t2)=26;
      elseif P_inipertrip(i2,j2,t2)>=28 && P_inipertrip(i2,j2,t2)<32
        P_inipertrip(i2,j2,t2)=30;
      elseif P_inipertrip(i2,j2,t2)>=32 && P_inipertrip(i2,j2,t2)<36
        P_inipertrip(i2,j2,t2)=34;
      elseif P_inipertrip(i2,j2,t2)>=36 && P_inipertrip(i2,j2,t2)<40
        P_inipertrip(i2,j2,t2)=38;
      elseif P_inipertrip(i2,j2,t2)>=40 && P_inipertrip(i2,j2,t2)<44
        P_inipertrip(i2,j2,t2)=42;
      elseif P_inipertrip(i2,j2,t2)>=44 && P_inipertrip(i2,j2,t2)<48
        P_inipertrip(i2,j2,t2)=46;
      elseif P_inipertrip(i2,j2,t2)>=48 && P_inipertrip(i2,j2,t2)<52
        P_inipertrip(i2,j2,t2)=50;
      elseif P_inipertrip(i2,j2,t2)>=52
        P_inipertrip(i2,j2,t2)=52;  
      end
      %}
    end
  end
end

% 122820-keep P_inipertrip ($/trip)

for t3=1:t_up 
     fid2 = fopen(['stdP24_Ppertripijt=',num2str(t3),'_061421','.txt'],'wt'); %creat & write Uijt to txt
     for i3=1:i_up
         for j3=1:j_up
             if j3 == j_up
                 fprintf(fid2,'%.2f\n',double(P_inipertrip(i3,j3,t3))); %0516-keep 4 decimals
             else
                 fprintf(fid2,'%.2f\t',double(P_inipertrip(i3,j3,t3))); %0516-keep 4 decimals
             end
         end    
     end    
     fclose(fid2);
end

%following is not working in CCR.
%{
for t21=1:40 %40 time periods
    fileName2=['Ini_Pijt=',num2str(t21),'_122820.txt']; %txt: Period0 to 39
    read2=load(fileName2); %load that txt
    P_inipertrip(:,:,t21)=read2; %assign txt to dijt
end
%}

% 111920-read Pijt
%{
for t2=1:40 %40 time periods
    fileName2=[num2str(Iter+1),'Iter_Pijt=',num2str(t2),'_112720.txt']; 
    read2=load(fileName2); %load that txt
    P(:,:,t2)=read2; %convert traveling time (min) to ($)
end
%}


% 0520-test whether z_op decreases as Pijt increases
%{
% for Pijt=1:1:1
%   z_op(Pijt) = zeros(Pijt);
% end
%}

%1-decision variable U111 to U(20,20,40),X111 to X(20,20,40),V1,1 to V(20,40)


% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
% X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
% %X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
% X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
% %X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
% U = sdpvar(i_up,j_up,t_up);
% %U = intvar(i_up,j_up,t_up); 
% U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
% %U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
% V = sdpvar(i_up,t_up); 
% %V = intvar(i_up,t_up); 
% V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
% %V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
% D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% %D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt

%1.1-parameters
% a=b=m=n=1,T=40
% a=1;
% b=1;
% m=1;
% n=1;
w_c=.88;%1.3;
w_w=1.08;%39;
w_p=.88;%1.3;
q=18.71;%58; %origin-2.62
k=10;%1; %origin-0.053
T=1; %total T is 40(look at obj func), this relates to t_up
T_abs=1; %temporarily set |T|=2 or 3 is related to game theory that |T|
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
% cijt=10;
pijt=1;%1.875; % 092920-this is $1.875/time period = $5/hrs
hit=5;%1;
alpha=1;

%initially price Pijt is very low, set to 1
%Pijt=1; %0520 - temporarily hide to see what happens as Pijt increases for oper_obj

% fid = fopen(['Pijt1to10001by1000PerStep_OperCostSplit.txt'],'wt'); %creat & write z_op to txt

%1.2-add operational level constraints

% Pijt=13;

% c1=0;

% CompU=zeros(16000,2);

% P_ij_SAV= zeros(i_up,j_up,t_up);

% for P_pa=1.3:0.1:1.4
%   P_pa 
% end


%100820-$6 to $24/TP = $17 to $68/hrs

% Pijt=13;

% c1=c1+1;

X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t

X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq

U = sdpvar(i_up,j_up,t_up);

U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 

V = sdpvar(i_up,t_up); 

V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0

D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt

E = sdpvar(i_up,j_up,t_up); %120520-empty trip Eijt

l1 = sdpvar(i_up,j_up,t_up); %102522-KKT Lagrange multipliers
l2 = sdpvar(i_up,j_up,t_up); %102522-KKT Lagrange multipliers
l3 = sdpvar(i_up,t_up); %102522-KKT Lagrange multipliers
l4 = sdpvar(t_up); %102522-KKT Lagrange multipliers

%{
multi=1.152;

for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if c(i,j,t)<=5
        P(i,j,t)=5.*multi;
      elseif c(i,j,t)>5 && c(i,j,t)<=10
        P(i,j,t)=10.*multi;
      elseif c(i,j,t)>10 && c(i,j,t)<=15
        P(i,j,t)=15.*multi;
      elseif c(i,j,t)>15 && c(i,j,t)<=20
        P(i,j,t)=20.*multi;
      elseif c(i,j,t)>20 && c(i,j,t)<=25
        P(i,j,t)=25.*multi;
      elseif c(i,j,t)>25 && c(i,j,t)<=30
        P(i,j,t)=30.*multi;
      elseif c(i,j,t)>30 && c(i,j,t)<=35
        P(i,j,t)=35.*multi;
      elseif c(i,j,t)>35 && c(i,j,t)<=40
        P(i,j,t)=40.*multi;
      elseif c(i,j,t)>40
        P(i,j,t)=45.*multi;
      end
    end
  end
end
%}

% 111820-avg cijt from t=1 to 40
%{
for t=1:t_up
  AvgC(t)=sum(sum(c(1:i_up,1:j_up,t)))/(i_up*j_up);
end
%}

% 111820-avg Pijt from t=1 to 40
%{
for t=1:t_up
  AvgP(t)=sum(sum(P(1:i_up,1:j_up,t)))/(i_up*j_up);
end
%}

% 100820-express D_ij_SAV
for i4=1:i_up
  for j4=1:j_up
    for t4=1:t_up
      D_ij_SAV(i4,j4,t4)=(.6.*c(i4,j4,t4).*d(i4,j4,t4).*w_c.*q-k.*w_w.*T_abs.*U(i4,j4,t4))./(k.*w_p.*P_inipertrip(i4,j4,t4)+c(i4,j4,t4).*w_c.*q);
%       Ddra(i,j,t) = D_ij_SAV(i,j,t)./d(i,j,t);
%       D_ij_SAV(i,j,t)=((k.*w_p.*P(i,j,t)+c(i,j,t).*w_c.*q).*T_abs.*U(i,j,t))./(c(i,j,t).*d(i,j,t).*w_c.*q-k.*w_w.*T_abs.*U(i,j,t));
%       D_ij_SAV(i,j,t)=((k.*w_p.*Pijt+cijt.*w_c.*q).*T_abs.*U(i,j,t))./(cijt.*d(i,j,t).*w_c.*q-k.*w_w.*T_abs.*U(i,j,t));
    end
  end
end

% 120520-express Eijt
for i5=1:i_up
  for j5=1:j_up
    for t5=1:t_up
      if t5 ==1
        E(i5,j5,t5)=X(i5,j5,t5)-D_ij_SAV(i5,j5,t5);
      else
        E(i5,j5,t5)=X(i5,j5,t5)-U(i5,j5,t5-1)-D_ij_SAV(i5,j5,t5);
      end
    end
  end
end

% %{

C = [];

% c26-092120overleaf
for i6=1:i_up %c26
   for j6=1:j_up
     for t6=1:t_up
       if (t6==1)
%          C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+Pijt))-X(i,j,1)];
         C = [C, U(i6,j6,t6) >= U_t0(i6,j6,t6)+D_ij_SAV(i6,j6,t6)-X(i6,j6,t6)]; %c26: t=1
       else       
         C = [C, U(i6,j6,t6) >= U(i6,j6,t6-1)+D_ij_SAV(i6,j6,t6)-X(i6,j6,t6)]; %c26: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i7=1:i_up %c27
    for j7=1:j_up
      for t7=1:t_up
        if (t7==1)
%           C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+Pijt))]; %c16: t=1
          C = [C, sum(X(i7,j7,t7:(t7+M))) >= U_t0(i7,j7,t7)+D_ij_SAV(i7,j7,t7)]; 
        else
          C = [C, sum(X(i7,j7,t7:(t7+M))) >= U(i7,j7,t7-1)+D_ij_SAV(i7,j7,t7)]; 
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i8=1:i_up %c28
  for t8=1:t_up
    if (t8==1)
      C = [C, V(i8,t8) == V_t0(i8,t8) + sum(alpha*X_t0(1:j_up,i8,t8)) - sum(X(i8,1:j_up,t8)) ]; 
    else
      C = [C, V(i8,t8) == V(i8,t8-1) + sum(sum(alpha*X(1:j_up,i8,1:(t8-1)))) - sum(X(i8,1:j_up,t8)) ]; %(remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c29-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i9=1:i_up %c30: Xijt
    for j9=1:j_up
        for t9=1:(t_up+M)
          C = [C, X(i9,j9,t9) >= 0];
        end
    end    
end    

for i10=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j10=1:j_up
        C = [C, X_t0(i10,j10,1) >= 0];
    end    
end    

for i11=1:i_up %c20: Uijt
    for j11=1:j_up
        for t11=1:t_up
            C = [C, U(i11,j11,t11) >= 0];
        end
    end    
end    

for i12=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j12=1:j_up
        C = [C, U_t0(i12,j12,1) >= 0];
    end    
end 

for i13=1:i_up %c20: Vit
    for t13=1:t_up
      C = [C, V(i13,t13) >= 0];
    end    
end    

for i14=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i14,1) >= 0];
end

% 092220-Uijt<=dijt

%for i15=1:i_up %c20: Uijt
%    for j15=1:j_up
%        for t15=1:t_up
%            C = [C, U(i15,j15,t15) <= d(i15,j15,t15)];
%        end
%    end    
%end    

% 102522-KKT-1st order condition for P-from Mathematica 102422
%{
for i21=1:i_up 
    for j21=1:j_up
        for t21=1:t_up
                      C = [C, (k.*(l1(i21,j21,t21)-l2(i21,j21,t21)+c(i21,j21,t21).*q)).*(.6.*c(i21,j21,t21).*d(i21,j21,t21).*q.*w_p-k.*T_abs.*U(i21,j21,t21).*w_w)./(k.*P_inipertrip(i21,j21,t21)+c(i21,j21,t21).*q)^2.*w_p == 0];
%             C = [C, (k.*(l1(i21,j21,t21)-l2(i21,j21,t21)+c(i21,j21,t21).*q)).*(.6.*c(i21,j21,t21).*d(i21,j21,t21).*q.*w_p-k.*T_abs.*U(i21,j21,t21).*w_w)./(k.*P_inipertrip(i21,j21,t21)+c(i21,j21,t21).*q)^2.*w_p == 0];
            %C = [C, X(i21,j21,t21)*U(i21,j21,t21)== 1000];
        end
    end    
end    
%}

%110122-KKT-1st order condition for X-from Mathematica 102422
%{
for i22=1:i_up 
    for j22=1:j_up
        for t22=1:t_up
%                       C = [C, l1(i22,j22,t22)-l2(i22,j22,t22) == 0];  
            C = [C, -c(i22,j22,t22)+l1(i22,j22,t22)-l2(i22,j22,t22)-l3(i22,j22)+l4(t22) == 0];  
        end
    end    
end  
%}

%110122-KKT-1st order condition for U-from Mathematica 102422
% %{
for i22=1:i_up 
    for j22=1:j_up
        for t22=1:t_up
            C = [C, -P_inipertrip(i22,j22,t22)-k.*l2(i22,j22,t22).*T_abs.*w_w./(k.*P_inipertrip(i22,j22,t22).*w_p+c(i22,j22,t22).*q.*w_p)-k.*P_inipertrip(i22,j22,t22).*T_abs.*w_w./(k.*P_inipertrip(i22,j22,t22).*w_p+c(i22,j22,t22).*q.*w_p)-l1(i22,j22,t22).*(-1-k.*T_abs.*w_w./(k.*P_inipertrip(i22,j22,t22).*w_p+c(i22,j22,t22).*q.*w_p))   == 0];  
        end
    end    
end  
%}

%110122-KKT-1st order condition for V-from Mathematica 102422
% %{
for i22=1:i_up 
    for j22=1:j_up
        for t22=1:t_up
            C = [C, -hit-l3(i22,j22)+l4(t22) == 0];  
        end
    end    
end  
%}

%110122-KKT-dual feasibility
for i9=1:i_up
    for j9=1:j_up
        for t9=1:t_up
          C = [C, l1(i9,j9,t9) >= 0];
          C = [C, l2(i9,j9,t9) >= 0];
        end
    end    
end    

% 110122-KKT-complementary slackness-lambda 1
%{
for i6=1:i_up 
   for j6=1:j_up
     for t6=1:t_up
       if (t6==1)
         C = [C,  l1(i6,j6,t6).*(U_t0(i6,j6,t6)-U(i6,j6,t6)+D_ij_SAV(i6,j6,t6)-X(i6,j6,t6))==0]; 
       else       
         C = [C,  l1(i6,j6,t6).*(U(i6,j6,t6-1)-U(i6,j6,t6)+D_ij_SAV(i6,j6,t6)-X(i6,j6,t6))==0]; 
       end
     end
   end    
end 
%}

% 110122-KKT-complementary slackness-lambda 2
%{
for i7=1:i_up 
    for j7=1:j_up
      for t7=1:t_up
        if (t7==1)
          C = [C,  l2(i7,j7,t7).*(U_t0(i7,j7,t7)-sum(X(i7,j7,t7:(t7+M)))+D_ij_SAV(i7,j7,t7))==0]; %
        else
          C = [C,  l2(i7,j7,t7).*(U(i7,j7,t7-1)-sum(X(i7,j7,t7:(t7+M)))+D_ij_SAV(i7,j7,t7))==0]; %
        end
      end  
    end     
end
%}


% 100820-Uijt<=Dijt
% for i=1:i_up 
%     for j=1:j_up
%         for t=1:t_up
%             C = [C, U(i,j,t) <= D_ij_SAV(i,j,t)];
%         end
%     end    
% end    

% 092920-Xijt<=dijt
% for i=1:i_up %c20: Uijt
%     for j=1:j_up
%         for t=1:t_up
%             C = [C, X(i,j,t) <= d(i,j,t)];
%         end
%     end    
% end    

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
%102522-by https://yalmip.github.io/command/sdpsettings/ YALMIP cannot solve KKT nonlinear calling CPLEX
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%ops = sdpsettings ('solver','cplex');
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(c(1:i_up,1:j_up,1:t_up).*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(c(1:i_up,1:j_up,1).*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));


% z_op(Pijt) = sum(sum(sum(c(i,j,t).*X(1:i_up,1:j_up,1:t_up+M)))) + sum(sum(c(i,j,t).*X_t0(1:i_up,1:j_up,1)))...
%   + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
%   + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
%   for i=1:i_up
%     for j=1:j_up
%       for t=1:t_up
%         D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+Pijt));
%       end
%     end
%   end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 101520
% result = optimize(C,z_op(Pijt),ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%   result = optimize(C,z);
if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

    fprintf('\nz_op(Pijt) value Iter=')
    value(z_op)
%     value(z_op(Pijt))
    
%     fprintf('\nPijt value: ')
%     value(P_pa)

    
%     fprintf(fid,'%.4f\t',double(z_op(Pijt))); %0523-save total cost;keep 4 decimals
%     fprintf(fid,'%.4f\t',double(sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1))))); 
%     fprintf(fid,'%.4f\t',double(sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))));
%     fprintf(fid,'%.4f\n',double(sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1))));
    
%   %     fprintf('\nD_ij_SAV value: ') %0515-as for now hide D_ij_SAV
%   %     value(D_ij_SAV) 
%   %  
%     fprintf('\nU value: ')
%     value(U)
    
    % 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid1 = fopen(['Pijt=',num2str(Pijt),'Uijt',num2str(U_t),'_092220','.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid1,'%.4f\n',double(U_count(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid1,'%.4f\t',double(U_count(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid1);
% end
    
%   %         
%       fprintf('\nX value: ') 
%       value(X)
%   %             
%   %     fprintf('\nV value: ')
%   %     value(V)
%   %     
%   %     fprintf('\nU_t0 value: ')
%   %     value(U_t0)
%   %     
%   %     fprintf('\nX_t0 value: ') 
%   %     value(X_t0)
%   %     

%   %     fprintf('\nV_t0 value: ')
%   %     value(V_t0)
% 
%   % else
%   % disp('Not solved!')
end


% 121920-keep Empty trip Eijt

for t16=1:t_up 
     fid1 = fopen([num2str(Iter),'stdP24_Eijt=',num2str(t16),'_061421','.txt'],'wt'); %creat & write Eijt
     for i16=1:i_up
         for j16=1:j_up
             if j16 == j_up
                 fprintf(fid1,'%.4f\n',double(E(i16,j16,t16))); %0516-keep 4 decimals
             else
                 fprintf(fid1,'%.4f\t',double(E(i16,j16,t16))); %0516-keep 4 decimals
             end
         end    
     end    
     fclose(fid1);
end

% 111820-keep Uijt
for t17=1:t_up 
     fid3 = fopen([num2str(Iter),'stdP24_Uijt=',num2str(t17),'_061421','.txt'],'wt'); %creat & write Uijt to txt
     for i17=1:i_up
         for j17=1:j_up
%            if U(i,j,U_t) <= 4.*D_ij_SAV(i,j,U_t)
             if j17 == j_up
                 fprintf(fid3,'%.4f\n',double(U(i17,j17,t17))); %0516-keep 4 decimals
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t)./d(i,j,U_t))); %0516-keep 4 decimals
             else
                 fprintf(fid3,'%.4f\t',double(U(i17,j17,t17))); %0516-keep 4 decimals
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t)./d(i,j,U_t))); %0516-keep 4 decimals
             end
%            end
         end    
     end    
     fclose(fid3);
end

% 112820-keep Xijt
for t18=1:t_up 
     fid4 = fopen([num2str(Iter),'stdP24_Xijt=',num2str(t18),'_061421','.txt'],'wt'); %creat & write Uijt to txt
     for i18=1:i_up
         for j18=1:j_up
%            if U(i,j,U_t) <= 4.*D_ij_SAV(i,j,U_t)
             if j18 == j_up
                 fprintf(fid4,'%.4f\n',double(X(i18,j18,t18))); %0516-keep 4 decimals
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t)./d(i,j,U_t))); %0516-keep 4 decimals
             else
                 fprintf(fid4,'%.4f\t',double(X(i18,j18,t18))); %0516-keep 4 decimals
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t)./d(i,j,U_t))); %0516-keep 4 decimals
             end
%            end
         end    
     end    
     fclose(fid4);
end

% 121920-keep Pijt
%{
for t=1:t_up 
     fid2 = fopen([num2str(Iter),'Iter_Pijt=',num2str(t),'_121920','.txt'],'wt'); %creat & write Uijt to txt
     for i=1:i_up
         for j=1:j_up
             if j == j_up
                 fprintf(fid2,'%.2f\n',double(P(i,j,t))); %0516-keep 4 decimals
             else
                 fprintf(fid2,'%.2f\t',double(P(i,j,t))); %0516-keep 4 decimals
             end
         end    
     end    
     fclose(fid2);
end
%}


% 121320-save D/d ratio
% %{
for t19=1:t_up 
     fid5 = fopen([num2str(Iter),'stdP24_Ddra=',num2str(t19),'_061421','.txt'],'wt'); %creat & write Uijt to txt
     for i19=1:i_up
         for j19=1:j_up
             if j19 == j_up
                 fprintf(fid5,'%.4f\n',double(D_ij_SAV(i19,j19,t19)./d(i19,j19,t19))); %0516-keep 4 decimals
             else
                 fprintf(fid5,'%.4f\t',double(D_ij_SAV(i19,j19,t19)./d(i19,j19,t19))); %0516-keep 4 decimals
             end
         end    
     end    
     fclose(fid5);
end


% 1.8-save Vit
fid6 = fopen([num2str(Iter),'stdP24_Vit_061421','.txt'],'wt'); %creat & write Vit to txt
for t20=1:t_up 
    for i20=1:i_up
      if i20==i_up
        fprintf(fid6,'%.4f\n',double(V(i20,t20))); %0517-keep 4 decimals
      else
        fprintf(fid6,'%.4f\t',double(V(i20,t20))); %0517-keep 4 decimals
      end
    end        
end
fclose(fid6);

% 111820-Save Avg Pijt for t=1...40
%{
fid1 = fopen([num2str(Iter),'Iter_AvgP=1to40_111820','.txt'],'wt'); %creat & write Uijt to txt
for t=1:t_up
   fprintf(fid1,'%.2f\t',double(AvgP(t))); %0516-keep 4 decimals
end    
fclose(fid1);
%}

% 111820-Save Avg Pijt for t=1...40
%{
fid2 = fopen([num2str(Iter),'Iter_AvgC=1to40_111820','.txt'],'wt'); %creat & write Uijt to txt
for t=1:t_up
   fprintf(fid2,'%.2f\t',double(AvgC(t))); %0516-keep 4 decimals
end    
fclose(fid2);
%}

% 0520-save total cost z_op(Pijt) in txt

% fclose(fid);

% 0517-following is save DijSAV in txt
%1.5-save D_ij_SAV to txt for plot use and find Pijt upper bound
% for t=1:t_up 
%     fid = fopen(['D_ij_SAVt',num2str(t),'Pijt=1Iter1.txt'],'wt'); %creat & write D_ij_SAV to txt
%     for i=1:i_up
%         for j=1:j_up
%             if j == j_up
%                 fprintf(fid,'%.4f\n',double(D_ij_SAV(i,j,t))); %0517-keep 4 decimals
%             else
%                 fprintf(fid,'%.4f\t',double(D_ij_SAV(i,j,t))); %0517-keep 4 decimals
%             end     
%         end    
%     end    
%     fclose(fid);
% end

% 0517-following is save X in txt, also include time window one
% 1.7-save Xijt
%{
for X_SaveTxt_t=1:(t_up+M)
  fid = fopen(['Xijt',num2str(X_SaveTxt_t),'Pijt=1Iter1.txt'],'wt'); %creat & write Xijt to txt
  for i=1:i_up
    for j=1:j_up
      if j==j_up
        fprintf(fid,'%.4f\n',double(X(i,j,X_SaveTxt_t))); %0517-keep 4 decimals
      else
        fprintf(fid,'%.4f\t',double(X(i,j,X_SaveTxt_t))); %0517-keep 4 decimals
      end     
    end    
   end    
   fclose(fid);
end
%}

% 0517-following is save V in txt
% 1.8-save Vit
% fid = fopen(['VitAllPijt=1Iter1.txt'],'wt'); %creat & write Vit to txt
% for V_SaveTxt_t=1:t_up 
%     for i=1:i_up
%       if i==i_up
%         fprintf(fid,'%.4f\n',double(V(i,V_SaveTxt_t))); %0517-keep 4 decimals
%       else
%         fprintf(fid,'%.4f\t',double(V(i,V_SaveTxt_t))); %0517-keep 4 decimals
%       end
%     end        
% end
% fclose(fid);

% 0519-temporarily hide Vi_t0
% 1.9-save Vi_t0 in txt
% fid = fopen(['Vi_t0Pijt=1.txt'],'wt'); %creat & write Vit to txt
% for i=1:i_up
%   fprintf(fid,'%.4f\t',double(V_t0(i,1))); %0517-keep 4 decimals
% end
% fclose(fid);


%} 

%0-end recording running time
toc